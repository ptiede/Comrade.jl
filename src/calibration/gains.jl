using LinearAlgebra
using SparseArrays
import Distributions
using Statistics
using PrettyTables

export GainCache, GainPrior, GainModel

"""
    $(TYPEDEF)

Abstract type that encompasses all RIME style corruptions.
"""
abstract type RIMEModel <: AbstractModel end

basemodel(m::RIMEModel) = m.model
flux(m::RIMEModel) = flux(basemodel(m))
radialextent(m::RIMEModel) = radialextent(basemodel(m))

function intensitymap(model::RIMEModel, fov::NTuple{2}, dim::Dims{2}, args...; kwargs...)
    return intensitymap(basemodel(model), fov, dim, args...; kwargs...)
end

"""
    $(TYPEDEF)

Internal type that holds the gain design matrices for visibility corruption.
See [`GainCache`](@ref GainCache).
"""
struct DesignMatrix{X,M<:AbstractMatrix{X},T,S} <: AbstractMatrix{X}
    matrix::M
    times::T
    stations::S
end

Base.getindex(m::DesignMatrix, i::Int) = getindex(m.matrix, i)
Base.size(m::DesignMatrix) = size(m.matrix)
Base.IndexStyle(::Type{<:DesignMatrix{X,M}}) where {X,M} = Base.IndexStyle(M)
Base.getindex(m::DesignMatrix, I::Vararg{Int,N}) where {N} = getindex(m.matrix, I...)
Base.setindex!(m::DesignMatrix, v, i::Int) = setindex!(m.matrix, v, i)
Base.setindex!(m::DesignMatrix, v, i::Vararg{Int, N}) where {N} = setindex!(m.matrix, v, i...)

Base.similar(m::DesignMatrix, ::Type{S}, dims::Dims) where {S} = DesignMatrix(similar(m.matrix, S, dims), m.times, m.stations)




function LinearAlgebra.mul!(y::AbstractArray, M::DesignMatrix, x::AbstractArray)
    LinearAlgebra.mul!(y, M.matrix, x)
end


struct GainPrior{S,T,G} <: Distributions.ContinuousMultivariateDistribution
    stations::S
    times::T
    dist::G
end

"""
    GainPrior(dists, st::ScanTable)

Creates a distribution for the gain priors for scan table `st`. The `dists` should be
a NamedTuple of `Distributions`, where each name corresponds to a telescope or station
in the scan table. See [`scantable`](@ref scantable). The resulting type if a subtype
of the `Distributions.AbstractDistribution` so the usual `Distributions` interface
should work.

# Example

For the 2017 observations of M87 a common GainPrior call is:
```julia-repl
julia> gdist = GainPrior((AA = LogNormal(0.0, 0.1),
                   AP = LogNormal(0.0, 0.1),
                   JC = LogNormal(0.0, 0.1),
                   SM = LogNormal(0.0, 0.1),
                   AZ = LogNormal(0.0, 0.1),
                   LM = LogNormal(0.0, 1.0),
                   PV = LogNormal(0.0, 0.1)
                ), st)

julia> x = rand(gdist)
julia> logdensityof(gdist, x)
```
"""
function GainPrior(dists, st::ScanTable)
    @argcheck Set(keys(dists)) == Set(stations(st))

    gtimes, gstat = gain_stations(st)

    gprior = Dists.product_distribution([getproperty(dists, g) for g in gstat])
    return GainPrior(gstat, gtimes, gprior)
end

#HypercubeTransform.bijector(d::GainPrior) = HypercubeTransform.asflat(d.dist)
HypercubeTransform.asflat(d::GainPrior) = asflat(d.dist)
HypercubeTransform.ascube(d::GainPrior) = ascube(d.dist)

Distributions.sampler(d::GainPrior) = Distributions.sampler(d.dist)
Base.length(d::GainPrior) = length(d.dist)
Base.eltype(d::GainPrior) = eltype(d.dist)
function Distributions._rand!(rng::AbstractRNG, d::GainPrior, x::AbstractVector)
    Distributions._rand!(rng, d.dist, x)
end

function Distributions._logpdf(d::GainPrior, x::AbstractArray)
    return Distributions._logpdf(d.dist, x)
end

struct HeirarchicalGainPrior{G,DM,DS,S,T}
    mean::DM
    std::DS
    stations::S
    times::T
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

HypercubeTransform.asflat(d::NamedDist{N}) where {N} = asflat(NamedTuple{N}(d.dists))

DensityInterface.DensityKind(::NamedDist) = DensityInterface.IsDensity()
DensityInterface.logdensityof(d::NamedDist, x) = Dists.logpdf(d, x)

DensityInterface.DensityKind(::HeirarchicalGainPrior) = DensityInterface.IsDensity()
DensityInterface.logdensityof(d::HeirarchicalGainPrior, x) = Dists.logpdf(d, x)

function _construct_gain_prior(means::NamedTuple{N}, stds::NamedTuple{N}, ::Type{G}, stations) where {N, G}
    gpr = NamedTuple{N}(G.(values(means), values(stds)))
    return Dists.product_distribution(getproperty.(Ref(gpr), stations))
end


function HeirarchicalGainPrior{G}(means, std, st::ScanTable) where {G}
    mnt = NamedDist(means)
    snt = NamedDist(std)
    gtimes, gstat = gain_stations(st)

    return HeirarchicalGainPrior{G, typeof(mnt), typeof(snt), typeof(gstat), typeof(gtimes)}(mnt, snt, gstat, gtimes)
end

function Dists.logpdf(d::HeirarchicalGainPrior{G}, x::NamedTuple) where {G}
    lm = Dists.logpdf(d.mean, x.mean)
    ls = Dists.logpdf(d.std, x.std)
    dg = _construct_gain_prior(x.mean, x.std, G, d.stations)
    lg = Dists.logpdf(dg, x.gains)
    return lg+ls+lm
end

function _unwrapped_logpdf(d::HeirarchicalGainPrior, x::Tuple)
    return Dists.logpdf(d, NamedTuple{(:mean, :std, :gains)}(x))
end


function Dists.rand(rng::AbstractRNG, d::HeirarchicalGainPrior{G}) where {G}
    m = rand(rng, d.mean)
    s = rand(rng, d.std)
    dg = _construct_gain_prior(m, s, G, d.stations)
    g = rand(rng, dg)
    return (mean=m, std=s, gains=g)
end

Base.length(d::HeirarchicalGainPrior) = length(d.mean) + length(d.std) + length(d.times)

function HypercubeTransform.asflat(d::HeirarchicalGainPrior{G}) where {G}
    m = rand(d.mean)
    s = rand(d.std)
    dg = _construct_gain_prior(m, s, G, d.stations)
    return TransformVariables.as((mean = asflat(d.mean), std = asflat(d.std), gains = asflat(dg)))
end

Statistics.mean(d::GainPrior) = mean(d.dist)
Statistics.var(d::GainPrior) = var(d.dist)
Distributions.entropy(d::GainPrior) = entropy(d.dist)
Statistics.cov(d::GainPrior) = cov(d.dist)

"""
    $(TYPEDEF)

# Fields
$(FIELDS)

# Notes
Internal type. This should be created using the [`GainCache(st::ScanTable)`](@ref GainCache) method.
"""
struct GainCache{D1<:DesignMatrix, D2<:DesignMatrix, T, S}
    """
    Gain design matrix for the first station
    """
    m1::D1
    """
    Gain design matrix for the second station
    """
    m2::D2
    """
    Set of times for each gain
    """
    times::T
    """
    Set of stations for each gain
    """
    stations::S
end

"""
    caltable(obs::EHTObservation, gains::AbstractVector)

Create a calibration table for the observations `obs` with `gains`. This returns
a [`CalTable`](@ref) object that satisfies the
`Tables.jl` interface. This table is very similar to the `DataFrames` interface.

# Example

```julia
ct = caltable(obs, gains)

# Access a particular station (here ALMA)
ct[:AA]
ct.AA

# Access a the first row
ct[1, :]
```
"""
function caltable(obs::EHTObservation, gains::AbstractVector)
    st = scantable(obs)
    gcache = GainCache(st)
    return caltable(gcache, gains)
end

"""
    caltable(g::GainCache, gains::AbstractVector)

Convert the GainCache `g` and recovered `gains` into a `CalTable` which satisfies the
`Tables.jl` interface. This table is very similar to the `DataFrames` interface.

# Example

```julia
ct = caltable(gcache, gains)

# Access a particular station (here ALMA)
ct[:AA]
ct.AA

# Access a the first row
ct[1, :]
```
"""
function caltable(g::GainCache, gains::AbstractVector)
    @argcheck length(g.times) == length(gains)

    stations = sort(unique(g.stations))
    times = unique(g.times)
    gmat = Matrix{Union{eltype(gains), Missing}}(missing, length(times), length(stations))

    alltimes = g.times
    allstations = g.stations
    # Create the lookup dict
    lookup = Dict(stations[i]=>i for i in eachindex(stations))
    for i in eachindex(gains)
        t = findfirst(t->(t==alltimes[i]), times)
        c = lookup[allstations[i]]
        gmat[t,c] = gains[i]
    end

    return CalTable(stations, lookup, times, gmat)
end


"""
    GainCache(st::ScanTable)

Creates a cache for the application of gain corruptions to the model visibilities.
This cache consists of the gain design matrices for each station and the set of times
and stations for each gain.

"""
function GainCache(st::ScanTable)
    gtime, gstat = gain_stations(st)
    m1, m2 = gain_design(st)
    return GainCache(m1, m2, gtime, gstat)
end

"""
    $(TYPEDEF)

A model that applies gain corruptions to a `Comrade` `model`.
This obeys the usual `Comrade` interface and can be evaluated using
`visibilities`.

# Fields
$(FIELDS)
"""
struct GainModel{C, G<:AbstractArray, M} <: RIMEModel
    """
    Cache for the application of gain. This can be constructed with
    [`GainCache`](@ref GainCache).
    """
    cache::C
    """
    Array of the specific gains that are to be applied to the visibilities.
    """
    gains::G
    """
    Base model that will be used to compute the uncorrupted visibilities.
    """
    model::M
end

"""
    caltable(g::GainModel)

Compute the gain calibration table from the [`GainModel`](@ref) `g`. This will
return a [`CalTable`](@ref) object, whose rows are different times,
and columns are different telescopes.
"""
function caltable(g::GainModel)
    return caltable(g.cache, g.gains)
end


function intensitymap!(img::IntensityMap, model::GainModel, p)
    return intensitymap!(img, model.model, p)
end

function _visibilities(model::GainModel, p)
    vis = _visibilities(model.model, p)
    return corrupt(vis, model.cache, model.gains)
end

function amplitudes(model::GainModel, p)
    amp = amplitudes(model.model, p)
    return abs.(corrupt(amp, model.cache, model.gains))
end

# Pass through since closure phases are independent of gains
function closure_phases(model::GainModel, args::Vararg{<:AbstractArray, N}) where {N}
    return closure_phases(model.model, args...)
end


# Pass through since log-closure amplitudes are independent of gains
function logclosure_amplitudes(model::GainModel, args::Vararg{<:AbstractArray, N}) where {N}
    return logclosure_amplitudes(model.model, args...)
end

"""
    corrupt(vis::AbstractArray, cache::GainCache, gains::AbstractArray)

Corrupt the visibilities `vis` with the gains `gains` using a `cache`.

This returns an array of corrupted visibilties. This is called internally
by the `GainModel` when producing the visibilties.
"""
function corrupt(vis::AbstractArray, cache::GainCache, gains::AbstractArray)
    g1 = cache.m1.matrix*gains
    g2 = cache.m2.matrix*gains
    return @. g1*vis*conj(g2)
end

# ChainRulesCore.@non_differentiable getproperty(cache::GainCache, s::Symbol)

# function ChainRulesCore.rrule(::typeof(corrupt), vis::AbstractArray, cache::GainCache, gains::AbstractArray)
#     g1 = cache.m1*gains
#     cg2 = conj.(cache.m2*gains)
#     viscor = @. g1*vis*cg2
#     function _corrupt_pullback(ΔV)
#         cΔV = conj.(ΔV)
#         Δf = NoTangent()
#         Δvis   = @thunk(cΔV.*g1.*cg2)
#         Δcache = NoTangent()

#         tmp1 = Diagonal(vis.*g1)*cache.m1
#         tmp2 = Diagonal(vis.*cg2)*cache.m2
#         Δgains = ΔV'*tmp1 + ΔV'*tmp2
#         return (Δf, Δvis, Δcache, Δgains)
#     end
#     return viscor, _corrupt_pullback
# end


# This is an internal function that computes the set of stations from a ScanTable
function gain_stations(st::ScanTable)
    gainstat = Symbol[]
    times = eltype(st.times)[]
    for i in 1:length(st)
        s = stations(st[i])
        append!(gainstat, s)
        append!(times, fill(st[i].time, length(s)))
    end
    return times, gainstat
end


# Construct a gain design matrices for each baseline station in a scan table
function gain_design(st::ScanTable)

    # Construct the indices that will tag where gains are
    rowInd1 = Int[]
    colInd1 = Int[]
    rowInd2 = Int[]
    colInd2 = Int[]
    times = st.obs[:time]
    bl = st.obs[:baseline]
    gaintime, gainstat = gain_stations(st)
    gts = collect(zip(gaintime, gainstat))
    for i in 1:length(times)
        t = times[i]
        s1, s2 = bl[i]
        # now get the columns that corresponds to the gain
        c1 = findall(x->((x[1]==t)&&(x[2]==s1)), gts)
        c2 = findall(x->((x[1]==t)&&(x[2]==s2)), gts)
        append!(colInd1, c1)
        append!(rowInd1, fill(i, length(c1)))
        append!(colInd2, c2)
        append!(rowInd2, fill(i, length(c2)))
    end
    z = fill(1.0, length(rowInd1))
    m1 = sparse(rowInd1, colInd1, z, length(times), length(gaintime))
    m2 = sparse(rowInd2, colInd2, z, length(times), length(gaintime))
    return DesignMatrix(m1, gaintime, gainstat), DesignMatrix(m2, gaintime, gainstat)
end
