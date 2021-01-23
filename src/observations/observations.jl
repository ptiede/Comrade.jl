abstract type AbstractInterferometryDatum{T} end

abstract type AbstractVisibilityDatum{T} <: AbstractInterferometryDatum{T} end
abstract type AbstractLinearPolDatum{S<:AbstractVisibilityDatum, T} <: AbstractInterferometryDatum{T} end
abstract type AbstractCrossPolDatum{S,T} <: AbstractInterferometryDatum{T} end

abstract type ClosureProducts{T} <: AbstractInterferometryDatum{T} end

abstract type Observation{T} end


Base.@kwdef struct EHTObservation{F,T<:AbstractVisibilityDatum{F}, N} <: Observation{F}
    data::StructArray{T,1}
    mjd::N
    ra::F
    dec::F
    source::Symbol
    timetype::Symbol = :UTC
end


Base.@kwdef struct EHTVisibilityDatum{T<:Number} <: AbstractVisibilityDatum{T}
    vis::Complex{T}
    error::T
    u::T
    v::T
    time::T
    frequency::T
    bandwidth::T
    baselines::NTuple{2,Symbol}
end


Base.@kwdef struct EHTVisibilityAmplitude{T<:Number} <: AbstractVisibilityDatum{T}
    amp::T
    error::T
    u::T
    v::T
    time::T
    frequency::T
    bandwidth::T
    baselines::NTuple{2,Symbol}
end


Base.@kwdef struct EHTClosurePhaseDatum{T<:Number} <: ClosureProducts{T}
    phase::T
    error::T
    u1::T
    v1::T
    u2::T
    v2::T
    time::T
    frequency::T
    bandwidth::T
    baselines::NTuple{3,Symbol}
end
