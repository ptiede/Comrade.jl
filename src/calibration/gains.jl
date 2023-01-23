
"""
    $(TYPEDEF)

Abstract type that encompasses all RIME style corruptions.
"""
abstract type RIMEModel <: AbstractModel end

basemodel(m::RIMEModel) = m.model
flux(m::RIMEModel) = flux(basemodel(m))
radialextent(m::RIMEModel) = radialextent(basemodel(m))

function intensitymap(model::RIMEModel, dims::AbstractDims)
    return intensitymap(basemodel(model), dims)
end

function intensity_point(model::RIMEModel, p)
    return intensity_point(basemodel(model), p)
end

_amplitudes(m::RIMEModel, u, v, t, f) = abs.(_visibilities(m, u, v, t, f))


"""
    $(TYPEDEF)

Internal type that holds the gain design matrices for visibility corruption.
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
