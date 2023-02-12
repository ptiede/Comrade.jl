export CalPrior, HeirarchicalCalPrior

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
function CalPrior(dists::NamedTuple, jcache::JonesCache, reference=:none)
    gstat = jcache.stations
    @argcheck Set(keys(dists)) == Set(gstat)

    if reference === :none
        gprior = make_gdist(dists, gstat, jcache)
    else
        gprior = make_reference_gdist(dists, gstat, jcache, reference)
    end
    return CalPrior{typeof(gprior), typeof(jcache)}(gprior, jcache)
end

"""
    CalPrior(dist0::NamedTuple, dist_transition::NamedTuple, jcache::RelativeJonesCache)

Constructs a calibration prior in two steps. The first two arguments have to be a named tuple
of distributions, where each name corresponds to a site. The first argument is gain prior
for the first time stamp. The second argument is the relative gain prior for each subsequent
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
function CalPrior(dist0::NamedTuple, distt::NamedTuple, jcache::RelativeJonesCache)
    sites = Tuple(unique(jcache.stations))
    @argcheck Set(keys(dist0)) == Set(sites)
    @argcheck Set(keys(distt)) == Set(sites)

    times = jcache.times
    stations = jcache.stations
    stimes = NamedTuple{sites}(map(x->times[findall(==(x), stations)], sites))
    cdist = map(zip(jcache.stations, jcache.times)) do (g, t)
        t > first(getproperty(stimes, g)) && return getproperty(distt, g)
        return getproperty(dist0, g)
    end

    gprior = Dists.product_distribution(cdist)
    return CalPrior{typeof(gprior), typeof(jcache)}(gprior, jcache)
end


function make_gdist(dists, gstat, jcache)
    gg = map(p->getproperty(dists, p), gstat)
    return Dists.product_distribution(gg)
end

function make_reference_gdist(dists, gstat, jcache, refprior)
    list = _makelist(dists, gstat, jcache, refprior)
    return Dists.product_distribution(list)
end


function _makelist(dists, gstat, jcache, refprior)
    sites = collect(Set(gstat))
    idx = 1
    times = jcache.times
    scantimes = unique(jcache.times)
    list = map(enumerate(scantimes)) do (i,t)
        # Select the reference station (we cycle through)
        inds = findall(==(t), times)
        ref, idxnew = _selectref(gstat[inds], sites, idx)
        idx = idxnew
        return gainlist_scan(gstat[inds], ref, dists, refprior)
    end
    return reduce(vcat, list)
end

function gainlist_scan(stations, ref, dists, refprior)
    return map(stations) do s
        if s == ref
            return refprior
        else
            return getproperty(dists, s)
        end
    end
end

function _selectref(stations, sites, idx)
    if sites[idx] ∈ stations
        return sites[idx], max((idx + 1)%(length(sites)+1), 1)
    else
        for i in (idx+1):(length(sites)+idx)
            idxnew = max(i%(length(sites)+1), 1)
            if sites[idxnew] ∈ stations
                return sites[idxnew], max((idxnew+1)%(length(sites)+1), 1)
            end
        end
        throw(AssertionError("No sites found"))
    end
end

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

struct HeirarchicalCalPrior{G,DM,DS,J}
    mean::DM
    std::DS
    jcache::J
end

struct NamedDist{Names, D}
    dists::D
end

function NamedDist(d::NamedTuple{N}) where {N}
    d = values(d)
    return NamedDist{N,typeof(d)}(d)
end

function Dists.logpdf(d::NamedDist{N}, x::NamedTuple{N}) where {N}
    vt = values(x)
    dists = d.dists
    sum(map((dist, acc) -> Dists.logpdf(dist, acc), dists, vt))
end

function Dists.logpdf(d::NamedDist{N}, x::NamedTuple{M}) where {N,M}
    xsub = select(x, N)
    return Dists.logpdf(d, xsub)
end

function Dists.rand(rng::AbstractRNG, d::NamedDist{N}) where {N}
    return NamedTuple{N}(map(x->rand(rng, x), d.dists))
end

function Dists.rand(rng::AbstractRNG, d::NamedDist{Na}, n::Dims) where {Na}
    map(CartesianIndices(n)) do I
        rand(rng, d)
    end
end

HypercubeTransform.asflat(d::NamedDist{N}) where {N} = asflat(NamedTuple{N}(d.dists))
HypercubeTransform.ascube(d::NamedDist{N}) where {N} = ascube(NamedTuple{N}(d.dists))

DensityInterface.DensityKind(::NamedDist) = DensityInterface.IsDensity()
DensityInterface.logdensityof(d::NamedDist, x) = Dists.logpdf(d, x)

DensityInterface.DensityKind(::HeirarchicalCalPrior) = DensityInterface.IsDensity()
DensityInterface.logdensityof(d::HeirarchicalCalPrior, x) = Dists.logpdf(d, x)

function _construct_gain_prior(means::NamedTuple{N}, stds::NamedTuple{N}, ::Type{G}, stations) where {N, G}
    gpr = NamedTuple{N}(G.(values(means), values(stds)))
    d = map(p->getproperty(gpr, p), stations)
    return Dists.product_distribution(d)
end


function HeirarchicalCalPrior{G}(means, std, jcache::JonesCache) where {G}
    mnt = NamedDist(means)
    snt = NamedDist(std)

    return HeirarchicalCalPrior{G, typeof(mnt), typeof(snt), typeof(jcache)}(mnt, snt, jcache)
end

function Dists.logpdf(d::HeirarchicalCalPrior{G}, x::NamedTuple) where {G}
    lm = Dists.logpdf(d.mean, x.mean)
    ls = Dists.logpdf(d.std, x.std)
    dg = _construct_gain_prior(x.mean, x.std, G, stations(d.jcache))
    lg = Dists.logpdf(dg, x.gains)
    return lg+ls+lm
end

function _unwrapped_logpdf(d::HeirarchicalCalPrior, x::Tuple)
    return Dists.logpdf(d, NamedTuple{(:mean, :std, :gains)}(x))
end


function Dists.rand(rng::AbstractRNG, d::HeirarchicalCalPrior{G}) where {G}
    m = rand(rng, d.mean)
    s = rand(rng, d.std)
    dg = _construct_gain_prior(m, s, G, d.stations)
    g = rand(rng, dg)
    return (mean=m, std=s, gains=g)
end

Base.length(d::HeirarchicalCalPrior) = length(d.mean) + length(d.std) + length(d.times)

function HypercubeTransform.asflat(d::HeirarchicalCalPrior{G}) where {G}
    m = rand(d.mean)
    s = rand(d.std)
    dg = _construct_gain_prior(m, s, G, d.stations)
    return TransformVariables.as((mean = asflat(d.mean), std = asflat(d.std), gains = asflat(dg)))
end

Statistics.mean(d::CalPrior) = mean(d.dist)
Statistics.var(d::CalPrior) = var(d.dist)
Distributions.entropy(d::CalPrior) = entropy(d.dist)
Statistics.cov(d::CalPrior) = cov(d.dist)
