struct ArrayPrior{D, A, R, C}
    default_dist::D
    override_dist::A
    refant::R
    centroid_station::C
end


"""
    ArrayPrior(default_dist; refant=NoReference(), kwargs...)

Construct a prior for an entire array of sites.

 - The `default_dist` is the default distribution for all sites. This may be an
[`IIDSitePrior`](@ref) or a time-correlated [`GaussMarkovSitePrior`](@ref).
 - Different priors for specified sites can be set using kwargs.
 - The `refant`  set the reference antennae to be used and is typically only done for priors that
correspond to gain phases.

# Example

```julia
p = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0, 0.1)); LM = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)) refant=SEFDReference())
```

means that every site has a normal prior with mean 0 and 0.1 std. dev. except LM which is mean
zero and unit std. dev. Finally the refant is using the [`SEFDReference`](@ref) scheme.
"""
function ArrayPrior(dist; refant = NoReference(), phase = nothing, kwargs...)
    phase === nothing || phase === false || throw(
        ArgumentError(
            "the `phase=true` cumulative reparameterization has been removed. Model phases " *
                "as unwrapped reals with a time-correlated prior instead, e.g. " *
                "`GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = ..., τ = ...); anchored = true)`."
        )
    )
    return ArrayPrior(dist, kwargs, refant, nothing)
end


# function site_priors(d::ArrayPrior, array)
#     return site_tuple(array, d.default_dist; d.override_dist...)
# end


struct ObservedArrayPrior{D, S} <: Distributions.ContinuousMultivariateDistribution
    dists::D
    sitemap::S
end
Base.eltype(d::ObservedArrayPrior) = eltype(d.dists)
Base.length(d::ObservedArrayPrior) = length(d.dists)
Dists.logpdf(d::ObservedArrayPrior, x::AbstractVector{<:Number}) = Dists.logpdf(d.dists, parent(x))
Dists._rand!(rng::Random.AbstractRNG, d::ObservedArrayPrior, x::AbstractArray{<:Real}) = SiteArray(Dists._rand!(rng, d.dists, x), d.sitemap)
# Flat (`TVFlat`) path: a TransformVariables instrument transform wrapping the raw
# TV node of the inner distribution, so it slots into the single flat TV tree.
function PT.transport_node(d::ObservedArrayPrior, space::PT.TVFlat)
    inner = PT.transport_node(d.dists, space)
    return InstrumentTransform(inner, d.sitemap)
end

# Std-space (`StdUniform`/`StdNormal`) path: a ProbabilityTransports node wrapping the
# PT node of the inner distribution.
function PT.transport_node(d::ObservedArrayPrior, space::PT.AbstractStdDist)
    inner = PT.transport_node(d.dists, space)
    return StdInstrumentTransform(inner, d.sitemap)
end

function build_sitemap(d::ArrayPrior, array)
    # construct the site by site prior
    sites_prior = site_tuple(array, d.default_dist; d.override_dist...)

    T = array[:Ti]
    F = array[:Fr]
    bls = array[:sites]

    # Per site, find the observed (time-segment, channel) pairs of its segmentation
    # grids by binning each of the site's data points into its containing cell
    # (O(n log nseg); segments within one grid are disjoint, so the binary search over
    # the sorted grid is valid). The `Segment`s are construction-transient: their
    # explicit half-open (lo, hi) edges are the keys here and go straight into the
    # lookup's axes.
    tkeys = Tuple{Float64, Float64}[]
    fkeys = Tuple{Float64, Float64}[]
    slist = Symbol[]
    for s in keys(sites_prior)
        sp = sites_prior[s]
        inds_s = findall(x -> ((x[1] == s) || (x[2] == s)), bls)
        ts = T[inds_s]
        fs = F[inds_s]
        # merged multifrequency configs can repeat and unsort scan rows: sort for the
        # binary search (duplicated segments dedup through the key set)
        tseg = sort!(collect(timestamps(segmentation(sp), array)))
        fseg = sort!(collect(freqchannels(freqseg(sp), array)))
        seen = Set{NTuple{4, Float64}}()
        for k in eachindex(ts, fs)
            jt = searchsortedfirst(tseg, ts[k])
            (jt ≤ lastindex(tseg) && ts[k] ∈ tseg[jt]) || continue
            jf = searchsortedfirst(fseg, fs[k])
            (jf ≤ lastindex(fseg) && fs[k] ∈ fseg[jf]) || continue
            tk = (float(tseg[jt].lo), float(tseg[jt].hi))
            fk = (float(fseg[jf].lo), float(fseg[jf].hi))
            key = (tk..., fk...)
            key ∈ seen && continue
            push!(seen, key)
            push!(tkeys, tk)
            push!(fkeys, fk)
            push!(slist, s)
        end
    end

    # unique axes sorted by the segment edges
    tuni = sort!(unique(tkeys))
    funi = sort!(unique(fkeys))
    saxis = sort!(unique(slist))
    tc0 = Int.(indexin(tkeys, tuni))
    fc0 = Int.(indexin(fkeys, funi))
    sc0 = Int.(indexin(slist, saxis))

    # order the flat storage (frequency, time, site) with sites fastest changing — the
    # layout the transports and reference-fixing assume
    ord = sortperm(collect(zip(fc0, tc0, sc0)))
    tcode = tc0[ord]
    fcode = fc0[ord]
    scode = sc0[ord]

    tlo = first.(tuni)
    thi = last.(tuni)
    flo = first.(funi)
    fhi = last.(funi)
    return _build_sitelookup(
        (tlo .+ thi) ./ 2, tlo, thi,
        (flo .+ fhi) ./ 2, flo, fhi,
        saxis, tcode, fcode, scode
    )
end

function ObservedArrayPrior(d::ArrayPrior, array::EHTArrayConfiguration)
    smap = build_sitemap(d, array)
    site_dists = site_tuple(array, d.default_dist; d.override_dist...)
    # Route on the site-prior interface (see independent.jl): all-product priors use the
    # fast vectorized product distribution; any chained prior (e.g. time-correlated
    # Gauss-Markov) switches every site to the per-chain machinery in gauss_markov.jl.
    all(sp -> structure(sp) isa ProductStructure, values(site_dists)) ||
        return build_chain_observed(d, site_dists, smap, array)
    dists = build_dist(site_dists, smap, array, d.refant, d.centroid_station)
    return ObservedArrayPrior(dists, smap)
end


# Flat (`TVFlat`) form: a TransformVariables transform that fixes some indices.
struct PartiallyFixedTransform{T, I, F} <: TV.AbstractTransform
    transform::T
    variate_index::I
    fixed_index::I
    fixed_values::F
end

TV.dimension(t::PartiallyFixedTransform) = TV.dimension(t.transform)

# Scatter the variate `y` and the fixed (reference) `fixed_values` back into the full
# parameter vector. Under Reactant this `setindex!`/broadcast lowers to a `stablehlo.scatter`,
# but the bounds check iterates the index vector with scalar `getindex`, which is disallowed
# while tracing. ComradeReactantExt overloads this for traced arrays to wrap it in
# `@allowscalar`. See https://github.com/EnzymeAD/Reactant.jl/issues/2960.
function fill_partially_fixed!(yfv, variate_index, fixed_index, y, fixed_values)
    yfv[variate_index] = y
    yfv[fixed_index] .= fixed_values
    return yfv
end

# Shared by the flat and Std transforms: allocate the full parameter vector and scatter
# the variate and fixed values into it.
function expand_partially_fixed(t, y)
    yfv = similar(y, length(t.variate_index) + length(t.fixed_index))
    return fill_partially_fixed!(yfv, t.variate_index, t.fixed_index, y, t.fixed_values)
end

function TV.transform_with(flag::TV.LogJacFlag, t::PartiallyFixedTransform, x, index)
    y, ℓ, index = TV.transform_with(flag, t.transform, x, index)
    return expand_partially_fixed(t, y), ℓ, index
end


function TV.inverse_at!(x::AbstractArray, index, t::PartiallyFixedTransform, y)
    return TV.inverse_at!(x, index, t.transform, y[t.variate_index])
end

TV.inverse_eltype(t::PartiallyFixedTransform, y::Type) = TV.inverse_eltype(t.transform, y)


# Std-space (`StdUniform`/`StdNormal`) form: a ProbabilityTransports node that fixes
# some indices. Mirrors the flat form but carries no Jacobian.
struct StdPartiallyFixedTransform{T, I, F} <: PT.AbstractTransport
    transform::T
    variate_index::I
    fixed_index::I
    fixed_values::F
end

PT.dimension(t::StdPartiallyFixedTransform) = PT.dimension(t.transform)

function PT.pfwd_step(t::StdPartiallyFixedTransform, x, index)
    y, index = PT.pfwd_step(t.transform, x, index)
    return expand_partially_fixed(t, y), index
end

function PT.pback_step!(y::AbstractVector, index, t::StdPartiallyFixedTransform, x)
    return PT.pback_step!(y, index, t.transform, @view(x[t.variate_index]))
end

PT.pback_eltype(t::StdPartiallyFixedTransform) =
    PT.pback_eltype(t.transform)


struct PartiallyConditionedDist{D <: Distributions.ContinuousMultivariateDistribution, I, F} <: Distributions.ContinuousMultivariateDistribution
    dist::D
    variate_index::I
    fixed_index::I
    fixed_values::F
end

Base.length(d::PartiallyConditionedDist) = length(d.variate_index) + length(d.fixed_index)
Base.eltype(d::PartiallyConditionedDist) = eltype(d.dist)


Distributions.sampler(d::PartiallyConditionedDist) = d

function Distributions.logpdf(d::PartiallyConditionedDist, x::AbstractVector)
    xv = @view parent(x)[d.variate_index]
    return Dists.logpdf(d.dist, xv)
end

function Distributions._rand!(rng::AbstractRNG, d::PartiallyConditionedDist, x::AbstractArray{<:Real})
    rand!(rng, d.dist, @view(x[d.variate_index]))
    # Now adjust the other indices
    x[d.fixed_index] .= d.fixed_values
    return x
end

function PT.transport_node(t::PartiallyConditionedDist, space::PT.TVFlat)
    return PartiallyFixedTransform(
        PT.transport_node(t.dist, space), t.variate_index, t.fixed_index, t.fixed_values
    )
end
function PT.transport_node(t::PartiallyConditionedDist, space::PT.AbstractStdDist)
    return StdPartiallyFixedTransform(
        PT.transport_node(t.dist, space), t.variate_index, t.fixed_index, t.fixed_values
    )
end

function build_dist(dists::NamedTuple, smap::SiteLookup, array, refants, centroid_station)
    ts = times(smap)
    ss = sites(smap)
    fixedinds, vals = reference_indices(array, smap, refants)

    if !(centroid_station isa Nothing)
        centstat = keys(centroid_station)
        vals = values(centroid_station)
        centroid1 = findfirst(==(centstat[1]), ss)
        centroid2 = findfirst(==(centstat[2]), ss)
        (centroid1 === nothing || centroid2 === nothing) && throw(ArgumentError("Centroid station not found in site list"))
        append!(fixedinds, [centroid1, centroid2])
        vals = append!(collect(vals), [vals[1], vals[2]])
    end

    variateinds = setdiff(eachindex(ts), fixedinds)
    dist = map(variateinds) do i
        getproperty(dists, ss[i]).dist
    end
    dist = Dists.product_distribution(dist)
    length(fixedinds) == 0 && return dist
    return PartiallyConditionedDist(dist, variateinds, fixedinds, vals)
end
