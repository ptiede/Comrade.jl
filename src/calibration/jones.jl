abstract type AbstractJones{T} <: AbstractMatrix{T} end



struct Gain{S, T<:Number} <: AbstractJones{T}
    g1::T
    g2::T
end

struct DTerm{S, T<:Number} <: AbstractJones{T}
    d1::T
    d2::T
end
