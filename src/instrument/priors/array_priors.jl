struct ArrayPrior{D, A, R, C}
    default_dist::D
    override_dist::A
    refant::R
    phase::Bool
    centroid_station::C
end


"""
    ArrayPrior(default_dist; refant=NoReference(), phase=false, kwargs...)

Construct a prior for an entire array of sites.

 - The `default_dist` is the default distribution for all sites. Currently only `IIDSitePrior` is supported.
 - Different priors for specified sites can be set using kwargs.
 - The `refant`  set the reference antennae to be used and is typically only done for priors that
correspond to gain phases.
 - The `phase` argument is a boolean that specifies if
the prior is for a `phase` or not. *The phase argument is experimental and we
recommend setting it to false currently.*

# Example

```julia
p = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0, 0.1)); LM = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)) refant=SEFDReference())
```

means that every site has a normal prior with mean 0 and 0.1 std. dev. except LM which is mean
zero and unit std. dev. Finally the refant is using the [`SEFDReference`](@ref) scheme.
"""
function ArrayPrior(dist; refant = NoReference(), phase = false, kwargs...)
    # if centroid_station isa Tuple{<:Symbol, <:Symbol}
    #     centroid_station = NamedTuple{centroid_station}((0.0, 0.0))
    # end
    return ArrayPrior(dist, kwargs, refant, phase, nothing)
end


# function site_priors(d::ArrayPrior, array)
#     return site_tuple(array, d.default_dist; d.override_dist...)
# end


struct ObservedArrayPrior{D, S} <: Distributions.ContinuousMultivariateDistribution
    dists::D
    sitemap::S
    phase::Bool
end
Base.eltype(d::ObservedArrayPrior) = eltype(d.dists)
Base.length(d::ObservedArrayPrior) = length(d.dists)
Dists.logpdf(d::ObservedArrayPrior, x::AbstractVector{<:Number}) = Dists.logpdf(d.dists, parent(x))
Dists._rand!(rng::Random.AbstractRNG, d::ObservedArrayPrior, x::AbstractArray{<:Real}) = SiteArray(Dists._rand!(rng, d.dists, x), d.sitemap)
# Flat (`StdFlat`) path: a TransformVariables instrument transform wrapping the raw
# TV node of the inner distribution, so it slots into the single flat TV tree.
function PT.transport_node(d::ObservedArrayPrior, space::PT.StdFlat)
    inner = PT.transport_node(d.dists, space)
    d.phase && return MarkovInstrumentTransform(inner, d.sitemap)
    return InstrumentTransform(inner, d.sitemap)
end

# Std-space (`StdUniform`/`StdNormal`) path: a ProbabilityTransports node wrapping the
# PT node of the inner distribution.
function PT.transport_node(d::ObservedArrayPrior, space)
    inner = PT.transport_node(d.dists, space)
    d.phase && return StdMarkovInstrumentTransform(inner, d.sitemap)
    return StdInstrumentTransform(inner, d.sitemap)
end

function build_sitemap(d::ArrayPrior, array)
    # construct the site by site prior
    sites_prior = site_tuple(array, d.default_dist; d.override_dist...)

    # Now we need all possible times to make sure we have all combinations
    T = array[:Ti]
    F = array[:Fr]


    # Ok to so this we are going to construct the schema first over sites.
    # At the end we may re-order depending on the schema ordering we want
    # to use.
    lists = map(keys(sites_prior)) do s
        seg = segmentation(sites_prior[s])
        # get all the indices where this site is present
        inds_s = findall(x -> ((x[1] == s)||x[2] == s), array[:sites])
        # Get all the times and frequencies for that site
        ts = T[inds_s]
        fs = F[inds_s]
        tfs = zip(ts, fs)
        # Now makes the acceptable time stamps given the segmentation
        tstamp = timestamps(seg, array)
        fchan = freqchannels(SpectralWindow(), array)
        # Now we find commonalities
        tf = Tuple{eltype(tstamp), eltype(fchan)}[]
        for f in fchan, t in tstamp
            if any(x -> (x[1] ∈ t && x[2] ∈ f), tfs) && ((!((t, f) ∈ tf)))
                push!(tf, (t, f))
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

function ObservedArrayPrior(d::ArrayPrior, array::EHTArrayConfiguration)
    smap = build_sitemap(d, array)
    site_dists = site_tuple(array, d.default_dist; d.override_dist...)
    dists = build_dist(site_dists, smap, array, d.refant, d.centroid_station)
    return ObservedArrayPrior(dists, smap, d.phase)
end


# Flat (`StdFlat`) form: a TransformVariables transform that fixes some indices.
struct PartiallyFixedTransform{T, I, F} <: TV.AbstractTransform
    transform::T
    variate_index::I
    fixed_index::I
    fixed_values::F
end

TV.dimension(t::PartiallyFixedTransform) = TV.dimension(t.transform)

function TV.transform_with(flag::TV.LogJacFlag, t::PartiallyFixedTransform, x, index)
    y, ℓ, index = TV.transform_with(flag, t.transform, x, index)
    yfv = similar(y, length(t.variate_index) + length(t.fixed_index))
    yfv[t.variate_index] = y
    yfv[t.fixed_index] .= t.fixed_values
    return yfv, ℓ, index
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

function PT.transport_step(t::StdPartiallyFixedTransform, x, index)
    y, index = PT.transport_step(t.transform, x, index)
    yfv = similar(y, length(t.variate_index) + length(t.fixed_index))
    yfv[t.variate_index] = y
    yfv[t.fixed_index] .= t.fixed_values
    return yfv, index
end

function PT.pullback_step!(y::AbstractVector, index, t::StdPartiallyFixedTransform, x)
    return PT.pullback_step!(y, index, t.transform, x[t.variate_index])
end

PT.pullback_eltype(t::StdPartiallyFixedTransform, ::Type{T}) where {T} =
    PT.pullback_eltype(t.transform, T)


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

function PT.transport_node(t::PartiallyConditionedDist, space::PT.StdFlat)
    return PartiallyFixedTransform(
        PT.transport_node(t.dist, space), t.variate_index, t.fixed_index, t.fixed_values
    )
end
function PT.transport_node(t::PartiallyConditionedDist, space)
    return StdPartiallyFixedTransform(
        PT.transport_node(t.dist, space), t.variate_index, t.fixed_index, t.fixed_values
    )
end

function build_dist(dists::NamedTuple, smap::SiteLookup, array, refants, centroid_station)
    ts = smap.times
    ss = smap.sites
    # fs = smap.frequencies
    fixedinds, vals = reference_indices(array, smap, refants)

    if !(centroid_station isa Nothing)
        centstat = keys(centroid_station)
        vals = values(centroid_station)
        centroid1 = findfirst(==(centstat[1]), ss)
        centroid2 = findfirst(==(centstat[2]), ss)
        centroid === nothing && throw(ArgumentError("Centroid station not found in site list"))
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
