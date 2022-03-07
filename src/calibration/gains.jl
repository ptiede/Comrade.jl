using LinearAlgebra
using SparseArrays
import Distributions
using Statistics

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

similar(m::DesignMatrix, ::Type{S}, dims::Dims) where {S} = DesignMatrix(similar(m.matrix, S, dims), m.times, m.stations)


const VisAmpDatum = Union{EHTVisibilityAmplitudeDatum, EHTVisibilityDatum}


function LinearAlgebra.mul!(y::AbstractArray, M::DesignMatrix, x::AbstractArray)
    LinearAlgebra.mul!(y, M.matrix, x)
end

struct GainPrior{S,T,G} <: Distributions.ContinuousMultivariateDistribution
    stations::S
    times::T
    dist::G
end

function GainPrior(dists, st::ScanTable)
    @argcheck Set(keys(dists)) == Set(stations(st))

    gtimes, gstat = gain_stations(st)

    gprior = Dists.product_distribution([getproperty(dists, g) for g in gstat])
    return GainPrior(gstat, gtimes, gprior)
end

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


struct GainCache{D1<:DesignMatrix, D2<:DesignMatrix, T, S}
    m1::D1
    m2::D2
    times::T
    stations::S
end

function GainCache(st::ScanTable)
    gtime, gstat = gain_stations(st)
    m1, m2 = gain_design(st)
    return GainCache(m1, m2, gtime, gstat)
end

"""
    `GainModel`
Construct the corruption model for the telescope gains
"""
struct GainModel{C, G, M} <: AbstractModel
    cache::C
    gains::G
    model::M
end

function intensitymap!(img::IntensityMap, model::GainModel)
    return intensitymap!(img, model.model)
end

function visibilities(model::GainModel, u::AbstractArray, v::AbstractArray)
    vis = visibilities(model.model, u, v)
    return corrupt(vis, model.cache, model.gains)
end

function corrupt(vis::AbstractArray, g::GainCache, gains::AbstractArray)
    g1 = g.m1*gains
    g2 = g.m2*gains
    return @. g1*vis*conj(g2)
end



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
