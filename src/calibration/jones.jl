abstract type AbstractJones{T} <: AbstractMatrix{T} end


mutable struct Gain{T<:Number} <: AbstractJones{T}
    g1::NamedTuple{:a,T}
    g2::NamedTuple{:b,T}
end

mutable struct DTerm{T<:Number} <: AbstractJones{T}
    d1::NamedTuple{:a, T}
    d2::NamedTuple{:b, T}
end
