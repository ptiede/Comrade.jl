export CalPrior, HierarchicalCalPrior

struct CalPrior{S,J<:AbstractJonesCache} <: Distributions.ContinuousMultivariateDistribution
    dists::S
    jcache::J
end

"""
    CalPrior(dists, cache::JonesCache, reference=:none)

Creates a distribution for the gain priors for gain cache `cache`. The `dists` should be
a NamedTuple of `Distributions`, where each name corresponds to a telescope or station
in the observation. The resulting type is a subtype
of the `Distributions.AbstractDistribution` so the usual `Distributions` interface
should work.

# Example

For the 2017 observations of M87 a common CalPrior call is:
```julia-repl
julia> gdist = CalPrior((AA = LogNormal(0.0, 0.1),
                   AP = LogNormal(0.0, 0.1),
                   JC = LogNormal(0.0, 0.1),
                   SM = LogNormal(0.0, 0.1),
                   AZ = LogNormal(0.0, 0.1),
                   LM = LogNormal(0.0, 1.0),
                   PV = LogNormal(0.0, 0.1)
                ), cache)

julia> x = rand(gdist)
julia> logdensityof(gdist, x)
```
"""
function CalPrior(dists::NamedTuple, jcache::AbstractJonesCache)
    gstat = jcache.schema.sites
    @argcheck issubset(Set(gstat), Set(keys(dists)))

    gprior = make_gdist(dists, gstat, jcache)
    return CalPrior{typeof(gprior), typeof(jcache)}(gprior, jcache)
end

"""
    CalPrior(dist0::NamedTuple, dist_transition::NamedTuple, jcache::SegmentedJonesCache)

Constructs a calibration prior in two steps. The first two arguments have to be a named tuple
of distributions, where each name corresponds to a site. The first argument is gain prior
for the first time stamp. The second argument is the segmented gain prior for each subsequent
time stamp. For instance, if we have
```julia
dist0 = (AA = Normal(0.0, 1.0), )
distt = (AA = Normal(0.0, 0.1), )
```
then the gain prior for first time stamp that AA obserserves will be `Normal(0.0, 1.0)`.
The next time stamp gain is the construted from
```
g2 = g1 + ϵ1
```
where `ϵ1 ~ Normal(0.0, 0.1) = distt.AA`, and `g1` is the gain from the first time stamp.
In other words `distt` is the uncorrelated transition probability when moving from timestamp
i to timestamp i+1. For the typical pre-calibrated dataset the gain prior on `distt` can be
tighter than the prior on `dist0`.
"""
function CalPrior(dist0::NamedTuple, distt::NamedTuple, jcache::JonesCache)
    sites = Tuple(unique(jcache.schema.sites))
    gstat = jcache.schema.sites
    @argcheck issubset(Set(gstat), Set(keys(dist0)))
    @argcheck issubset(Set(gstat), Set(keys(distt)))

    times = jcache.schema.times
    stations = jcache.schema.sites
    stimes = NamedTuple{sites}(map(x->times[findall(==(x), stations)], sites))
    cdist = map(zip(jcache.schema.sites, jcache.schema.times)) do (g, t)
        ((t > first(getproperty(stimes, g))) && jcache.seg[g] isa ScanSeg{true} ) && return getproperty(distt, g)
        return getproperty(dist0, g)
    end

    gprior = Dists.product_distribution(cdist)
    return CalPrior{typeof(gprior), typeof(jcache)}(gprior, jcache)
end


function make_gdist(dists, gstat, jcache)
    gg = map(p->getproperty(dists, p), gstat)
    return Dists.product_distribution(gg)
end

# function make_reference_gdist(dists, gstat, jcache, refprior)
#     list = _makelist(dists, gstat, jcache, refprior)
#     return Dists.product_distribution(list)
# end


# function _makelist(dists, gstat, jcache, refprior)
#     sites = collect(Set(gstat))
#     idx = 1
#     times = jcache.schema.times
#     scantimes = unique(jcache.schema.times)
#     list = map(enumerate(scantimes)) do (i,t)
#         # Select the reference station (we cycle through)
#         inds = findall(==(t), times)
#         ref, idxnew = _selectref(gstat[inds], sites, idx)
#         idx = idxnew
#         return gainlist_scan(gstat[inds], ref, dists, refprior)
#     end
#     return reduce(vcat, list)
# end

# function gainlist_scan(stations, ref, dists, refprior)
#     return map(stations) do s
#         if s == ref
#             return refprior
#         else
#             return getproperty(dists, s)
#         end
#     end
# end

# function _selectref(stations, sites, idx)
#     if sites[idx] ∈ stations
#         return sites[idx], max((idx + 1)%(length(sites)+1), 1)
#     else
#         for i in (idx+1):(length(sites)+idx)
#             idxnew = max(i%(length(sites)+1), 1)
#             if sites[idxnew] ∈ stations
#                 return sites[idxnew], max((idxnew+1)%(length(sites)+1), 1)
#             end
#         end
#         throw(AssertionError("No sites found"))
#     end
# end

#HypercubeTransform.bijector(d::CalPrior) = HypercubeTransform.asflat(d.dist)
HypercubeTransform.asflat(d::CalPrior) = asflat(d.dists)
HypercubeTransform.ascube(d::CalPrior) = ascube(d.dists)

Distributions.sampler(d::CalPrior) = Distributions.sampler(d.dists)
Base.length(d::CalPrior) = length(d.dists)
Base.eltype(d::CalPrior) = eltype(d.dists)

function Distributions._rand!(rng::AbstractRNG, d::CalPrior, x::AbstractVector)
    Distributions._rand!(rng, d.dists, x)
end

function Distributions._logpdf(d::CalPrior, x::AbstractArray)
    return Distributions._logpdf(d.dists, x)
end

# struct HierarchicalCalPrior{G,DM,DS,J}
#     mean::DM
#     std::DS
#     jcache::J
# end


# DensityInterface.DensityKind(::HierarchicalCalPrior) = DensityInterface.IsDensity()
# DensityInterface.logdensityof(d::HierarchicalCalPrior, x) = Dists.logpdf(d, x)

# function _construct_gain_prior(means::NamedTuple{N}, stds::NamedTuple{N}, ::Type{G}, stations) where {N, G}
#     gpr = NamedTuple{N}(map(G, values(means), values(stds)))
#     d = map(p->getproperty(gpr, p), stations)
#     return Dists.product_distribution(d)
# end


# function HierarchicalCalPrior{G}(means, std, jcache::JonesCache) where {G}
#     return HierarchicalCalPrior{G, typeof(means), typeof(std), typeof(jcache)}(means, std, jcache)
# end

# function Dists.logpdf(d::HierarchicalCalPrior{G}, x::NamedTuple) where {G}
#     lm = Dists.logpdf(d.mean, x.mean)
#     ls = Dists.logpdf(d.std, x.std)
#     dg = _construct_gain_prior(x.mean, x.std, G, stations(d.jcache))
#     lg = Dists.logpdf(dg, x.gains)
#     return lg+ls+lm
# end

# function _unwrapped_logpdf(d::HierarchicalCalPrior, x::Tuple)
#     return Dists.logpdf(d, NamedTuple{(:mean, :std, :gains)}(x))
# end


# function Dists.rand(rng::AbstractRNG, d::HierarchicalCalPrior{G}) where {G}
#     m = map(x->rand(rng, x), d.mean)
#     s = map(x->rand(rng, x), d.std)
#     dg = _construct_gain_prior(m, s, G, d.jcache.schema.sites)
#     g = rand(rng, dg)
#     return (mean=m, std=s, gains=g)
# end

# Base.length(d::HierarchicalCalPrior) = length(d.mean) + length(d.std) + length(d.times)

# function HypercubeTransform.asflat(d::HierarchicalCalPrior{G}) where {G}
#     m = rand(d.mean)
#     s = rand(d.std)
#     dg = _construct_gain_prior(m, s, G, d.jcache.schema.sites)
#     return TransformVariables.as((mean = asflat(d.mean), std = asflat(d.std), gains = asflat(dg)))
# end

Statistics.mean(d::CalPrior) = mean(d.dists)
Statistics.var(d::CalPrior) = var(d.dists)
Statistics.cov(d::CalPrior) = cov(d.dists)
