export SitesPrior, ObservedSitesPrior

abstract type AbstractSitePrior <: Distributions.ContinuousMultivariateDistribution end

abstract type AbstractInstrumentPrior <: Distributions.ContinuousMultivariateDistribution end

include("segmentation.jl")
include("independent.jl")
include("refant.jl")
include("instrument.jl")

struct SitesPrior{T, D, A, R}
    timecorr::T
    default_dist::D
    override_dist::A
    refant::R
end

function SitesPrior(corr, dist; refant=NoReference(), kwargs...)
    return SitesPrior(corr, dist, kwargs, refant)
end


struct ObservedSitesPrior{D, S} <: Distributions.ContinuousMultivariateDistribution
    dists::D
    sitemap::S
end
Base.eltype(d::ObservedSitesPrior) = eltype(d.dists)
Base.length(d::ObservedSitesPrior) = length(d.dists)
Dists._logpdf(d::ObservedSitesPrior, x::AbstractArray{<:Real}) = Dists._logpdf(d.dists, x)
Dists._rand!(rng::Random.AbstractRNG, d::ObservedSitesPrior, x::AbstractArray{<:Real}) = SiteArray(Dists._rand!(rng, d.dists, x), d.sitemap)
HypercubeTransform.asflat(d::ObservedSitesPrior) = asflat(d.dists)
HypercubeTransform.ascube(d::ObservedSitesPrior) = ascube(d.dists)

function build_sitemap(segmentation, array)
    sts = sites(array)
    ts  = timestamps(segmentation, array)
    fs  = unique(array[:F])

    T  = array[:T]
    bl = array[:sites]
    F  = array[:F]

    times = eltype(T)[]
    freqs = eltype(F)[]
    slist = Symbol[]

    # This is frequency, time, station ordering
    for f in fs, t in ts, s in sts
        if any(i->((T[i]∈t)&&s∈bl[i]&&f==F[i]), eachindex(T))
            push!(times, t.t0)
            push!(slist, s)
            push!(freqs, f)
        end
    end
    return SiteLookup(slist, times, freqs)
end


function ObservedSitesPrior(d::SitesPrior, array::EHTArrayConfiguration)
    smap = build_sitemap(d.timecorr, array)
    site_dists = sites_tuple(array, d.default_dist; d.override_dist...)
    dists = build_dist(site_dists, smap, array, d.refant)
    return ObservedSitesPrior(dists, smap)
end

function build_dist(dists::NamedTuple, smap::SiteLookup, array, refants)
    ts = smap.times
    ss = smap.sites
    # fs = smap.frequencies
    fixedinds, vals = reference_indices(array, smap, refants)

    variateinds = setdiff(eachindex(ts), fixedinds)
    dist = map(variateinds) do i
        getproperty(dists, ss[i])
    end
    dist = Dists.product_distribution(dist)
    length(fixedinds) == 0 && return dist
    return PartiallyConditionedDist(dist, variateinds, fixedinds, vals)
end
