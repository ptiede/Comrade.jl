using LinearAlgebra
using SparseArrays
import Distributions
using Statistics

export GainCache, GainPrior, GainModel

"""
    RIMEModel

Abstract type that encompasses all RIME style corruptions.
"""
abstract type RIMEModel <: AbstractModel end

basemodel(m::RIMEModel) = m.model
flux(m::RIMEModel) = flux(basemodel(m))
radialextent(m::RIMEModel) = radialextent(basemodel(m))

function intensitymap(model::RIMEModel, fovx::Number, fovy::Number, nx::Int, ny::Int, args...; kwargs...)
    return intensitymap(basemodel(model), fovx, fovy, nx, ny, args...; kwargs...)
end

"""
    DesignMatrix

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
Base.getindex(m::DesignMatrix, I::Vararg{Int,N}) where {N} = getindex(m.matrix, I)
Base.setindex!(m::DesignMatrix, v, i::Int) = setindex!(m.matrix, v, i)
Base.setindex!(m::DesignMatrix, v, i::Vararg{Int, N}) where {N} = setindex!(m.matrix, v, i)

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
    GainPrior(dists, st)

Creates a distribution for the gain priors for scan table `st`. The `dists` should be
a NamedTuple of `Distributions`, where each name corresponds to a telescope or station
in the scan table. See [`scantable`](@ref scantable). The resulting type if a subtype
of the `Distributions.AbstractDistribution` so the usual `Distributions` interface
should work.

# Example

For the 2017 observations of M87 a common GainPrior call is:
```julia
gdist = GainPrior((AA = LogNormal(0.0, 0.1),
                   AP = LogNormal(0.0, 0.1),
                   JC = LogNormal(0.0, 0.1),
                   SM = LogNormal(0.0, 0.1),
                   AZ = LogNormal(0.0, 0.1),
                   LM = LogNormal(0.0, 1.0),
                   PV = LogNormal(0.0, 0.1)
                ), st)

x = rand(gdist)
logdensityof(gdist, x)
```
"""
function GainPrior(dists, st::ScanTable)
    @argcheck Set(keys(dists)) == Set(stations(st))

    gtimes, gstat = gain_stations(st)

    gprior = Dists.product_distribution([getproperty(dists, g) for g in gstat])
    return GainPrior(gstat, gtimes, gprior)
end

HypercubeTransform.bijector(d::GainPrior) = HypercubeTransform.bijector(d.dist)
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
    $(SIGNATURES)

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


function intensitymap!(img::IntensityMap, model::GainModel, args...)
    return intensitymap!(img, model.model, args...)
end

function visibilities(model::GainModel, u::AbstractArray, v::AbstractArray)
    vis = visibilities(model.model, u, v)
    return corrupt(vis, model.cache, model.gains)
end

function amplitudes(model::GainModel, u::AbstractArray, v::AbstractArray)
    amp = amplitudes(model.model, u, v)
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
    $(SIGNATURES)

Corrupt the visibilities `vis` with the gains `gains` using a `cache`.

This returns an array of corrupted visibilties. This is called internally
by the `GainModel` when producing the visibilties.
"""
function corrupt(vis::AbstractArray, cache::GainCache, gains::AbstractArray)
    g1 = cache.m1*gains
    g2 = cache.m2*gains
    return @. g1*vis*conj(g2)
end


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
    z = fill(1, length(rowInd1))
    m1 = sparse(rowInd1, colInd1, z, length(times), length(gaintime))
    m2 = sparse(rowInd2, colInd2, z, length(times), length(gaintime))
    return DesignMatrix(m1, gaintime, gainstat), DesignMatrix(m2, gaintime, gainstat)
end
