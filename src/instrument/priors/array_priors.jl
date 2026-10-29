struct ArrayPrior{D, A, R, C, G, S}
    default_dist::D
    override_dist::A
    refant::R
    centroid_station::C
    gauge::G
    # the stations the parameter exists for; `nothing` for every station
    stations::S
end


"""
    ArrayPrior(default_dist; refant=NoReference(), gauge=nothing, stations=nothing, kwargs...)

Construct a prior for an entire array of sites.

 - The `default_dist` is the default distribution for all sites. This may be an
[`IIDSitePrior`](@ref) or a time-correlated [`GaussMarkovSitePrior`](@ref).
 - Different priors for specified sites can be set using kwargs named after the stations.
 - The `refant`  set the reference antennae to be used and is typically only done for priors that
correspond to gain phases.
 - `stations`, a tuple of station names, restricts the parameter to those stations, for a
term only some stations have (see [`PerStation`](@ref)); stations without data are skipped.
`nothing` (the default) is every station.
 - The `gauge` argument declares that this parameter is a term of a quantity whose absolute
level the data cannot see. `gauge = StationPhase()` marks the parameter as one summand of the
absolute station phase: the visibility phase on a baseline sees only differences of station
phases, so the sum of all `StationPhase()` parameters at a site is observable only up to one
constant per (time, frequency) stamp. When several such parameters are added together (a
track-level offset plus a per-integration residual, say) the reference pins placed on each
one separately can leave the *sum* under-determined. Declaring the terms lets the
[`InstrumentModel`](@ref) find the leftover flat directions with
[`Comrade.gauge_pins`](@ref) and either report them or pin them, depending on its
`gaugefix` setting. `gauge = nothing` (the default) opts out of that analysis.

# Example

```julia
p = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0, 0.1)); LM = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)), refant = SensitiveReference(0.0))
```

means that every site has a normal prior with mean 0 and 0.1 std. dev. except LM which is mean
zero and unit std. dev. Finally the refant is using the [`SensitiveReference`](@ref) scheme.
"""
function ArrayPrior(dist; refant = NoReference(), gauge = nothing, stations = nothing, kwargs...)
    # if centroid_station isa Tuple{<:Symbol, <:Symbol}
    #     centroid_station = NamedTuple{centroid_station}((0.0, 0.0))
    # end
    gauge isa Union{Nothing, StationPhase} ||
        throw(ArgumentError("ArrayPrior gauge must be nothing or StationPhase(), got $(repr(gauge))"))
    st = isnothing(stations) ? nothing : Tuple(Symbol.(stations))
    return ArrayPrior(dist, kwargs, refant, nothing, gauge, st)
end

# The site priors of `d` for the stations with measured values in `t` it covers.
function _site_priors(d::ArrayPrior, t::VisTable)
    s = _measuredpoints(t).sites
    present = sort(unique(vcat(first.(s), last.(s))))
    st = isnothing(d.stations) ? present : filter(in(present), collect(d.stations))
    isempty(st) && throw(ArgumentError("none of the stations $(d.stations) of an instrument term has data"))
    return site_tuple(st, d.default_dist; d.override_dist...)
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
    return InstrumentTransform(PT.transport_node(d.dists, space), d.sitemap)
end

# Std-space (`StdUniform`/`StdNormal`) path: a ProbabilityTransports node wrapping the
# PT node of the inner distribution.
function PT.transport_node(d::ObservedArrayPrior, space::PT.AbstractStdDist)
    return StdInstrumentTransform(PT.transport_node(d.dists, space), d.sitemap)
end

function build_sitemap(d::ArrayPrior, t::VisTable)
    # construct the site by site prior
    sites_prior = _site_priors(d, t)

    # Now we need all possible times to make sure we have all combinations
    points = _measuredpoints(t)
    T = points.Ti
    F = points.Fr
    fchan = [bandchannel(t)]


    # Ok to so this we are going to construct the schema first over sites.
    # At the end we may re-order depending on the schema ordering we want
    # to use.
    lists = map(keys(sites_prior)) do s
        seg = segmentation(sites_prior[s])
        # get all the indices where this site is present
        inds_s = findall(x -> ((x[1] == s)||x[2] == s), points.sites)
        # Get all the times and frequencies for that site
        ts = T[inds_s]
        fs = F[inds_s]
        tfs = zip(ts, fs)
        # Now makes the acceptable time stamps given the segmentation
        tstamp = timestamps(seg, t)
        # Now we find commonalities
        tf = Tuple{eltype(tstamp), eltype(fchan)}[]
        for f in fchan, τ in tstamp
            if any(x -> (x[1] ∈ τ && x[2] ∈ f), tfs) && ((!((τ, f) ∈ tf)))
                push!(tf, (τ, f))
            end
        end
        return first.(tf), last.(tf)
    end
    tlists = first.(lists)
    flists = last.(lists)
    # construct the schema
    slist = mapreduce((t, s) -> fill(s, length(t)), vcat, tlists, keys(sites_prior))
    tlist = reduce(vcat, tlists)
    flist = reduce(vcat, flists)


    tlistre = similar(tlist)
    slistre = similar(slist)
    flistre = similar(flist)
    # Now rearrange so we have frquency, time, site ordering (sites are the fastest changing)
    tuni = sort(unique((tlist)))
    funi = sort(unique((flist)))
    ind0 = 1
    for f in funi, t in tuni
        ind = (f .== flist) .& (t .== tlist)
        vtlist = @view tlist[ind]
        vslist = @view slist[ind]
        vflist = @view flist[ind]
        tlistre[ind0:(ind0 + length(vtlist) - 1)] .= vtlist
        slistre[ind0:(ind0 + length(vtlist) - 1)] .= vslist
        flistre[ind0:(ind0 + length(vtlist) - 1)] .= vflist
        ind0 += length(vtlist)
    end
    return SiteLookup(tlistre, flistre, slistre)
end

function ObservedArrayPrior(d::ArrayPrior, t::VisTable)
    smap = build_sitemap(d, t)
    site_dists = _site_priors(d, t)
    # Time-correlated site priors need per-chain machinery rather than a product
    # distribution; see gauss_markov.jl.
    any(needs_chain_machinery, values(site_dists)) &&
        return build_markov_observed(d, site_dists, smap, t)
    dists = build_dist(site_dists, smap, t, d.refant, d.centroid_station)
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

# Scatter `vals` into `y` at (static, host-side) `inds`. Under Reactant the broadcast
# lowers to a single vectorized `stablehlo.scatter` — the fast, "raised" path — instead of
# a serialized element-by-element traced loop; ComradeReactantExt overloads this for
# traced arrays (the bounds check scalar-iterates the index vector, and the
# single-element scatter hits a Reactant bug — see EnzymeAD/Reactant.jl#2960).
function scatter_values!(y, inds, vals)
    y[inds] .= vals
    return y
end

# Scatter the variate `y` and the fixed (reference) `fixed_values` back into the full
# parameter vector. Both scatters go through `scatter_values!`, which carries the traced
# workarounds (ComradeReactantExt overloads it — see EnzymeAD/Reactant.jl#2960), so this
# generic method serves both backends.
function fill_partially_fixed!(yfv, variate_index, fixed_index, y, fixed_values)
    scatter_values!(yfv, variate_index, y)
    scatter_values!(yfv, fixed_index, fixed_values)
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

function build_dist(dists::NamedTuple, smap::SiteLookup, t, refants, centroid_station)
    ts = smap.times
    ss = smap.sites
    # fs = smap.frequencies
    fixedinds, vals = reference_indices(t, smap, refants)
    if any(!isnothing ∘ initprior, values(dists))
        fixedinds = convert(Vector{Int}, fixedinds)
        vals = vals === nothing ? Float64[] : convert(Vector{Float64}, vals)
        # Init pins: a site prior with `init = FixedInit(v)` fixes each site's first
        # stamp to `v` — see `IIDSitePrior`. A stamp the
        # referencing scheme already fixes stays the reference's, and the values must
        # agree: a silent disagreement would quietly change the reference gauge.
        for s in keys(dists)
            ini = initprior(getproperty(dists, s))
            ini === nothing && continue
            for f in unique(smap.frequencies)
                sinds = findall(
                    i -> ss[i] == s && smap.frequencies[i] == f, eachindex(ss)
                )
                isempty(sinds) && continue
                i1 = sinds[argmin(view(ts, sinds))]
                j = findfirst(==(i1), fixedinds)
                if j === nothing
                    push!(fixedinds, i1)
                    push!(vals, ini.value)
                elseif vals[j] != ini.value
                    throw(
                        ArgumentError(
                            "init = FixedInit($(ini.value)) for site $s conflicts " *
                                "with the reference value $(vals[j]) already fixing " *
                                "its first stamp"
                        )
                    )
                end
            end
        end
    end

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
