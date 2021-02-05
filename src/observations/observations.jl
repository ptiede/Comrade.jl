abstract type AbstractInterferometryDatum{T} end

abstract type AbstractVisibilityDatum{T} <: AbstractInterferometryDatum{T} end
abstract type AbstractLinearPolDatum{S<:AbstractVisibilityDatum, T} <: AbstractInterferometryDatum{T} end
abstract type AbstractCrossPolDatum{S,T} <: AbstractInterferometryDatum{T} end

abstract type ClosureProducts{T} <: AbstractInterferometryDatum{T} end

abstract type Observation{T} end

using DelimitedFiles
using AstroTime: UTCEpoch, modified_julian

Base.@kwdef struct EHTObservation{F,T<:AbstractVisibilityDatum{F}, N} <: Observation{F}
    data::StructArray{T}
    mjd::N
    ra::F
    dec::F
    source::Symbol
    timetype::Symbol = :UTC
end

getdata(obs::Observation, s::Symbol) = getproperty(obs.data, s)
nsamples(obs) = length(obs.data)
Base.getindex(obs::Observation, i) = getindex(obs.data, i)


Base.@kwdef struct EHTVisibilityDatum{T<:Number} <: AbstractVisibilityDatum{T}
    visr::T
    visi::T
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

function load_tpy(file)
    data = readdlm(file, skipstart=1)
    bs1 = Symbol.(getindex.(data[:,6], Ref(1:2)))
    bs2 = Symbol.(getindex.(data[:,6], Ref(3:4)))
    baselines = tuple.(bs1, bs2)
    edata = StructArray{EHTVisibilityDatum{Float64}}(
                visr=float.(data[:,9]),
                visi=float.(data[:,11]),
                error=float.(data[:,10]),
                u=float.(data[:,7])*1e6,
                v=float.(data[:,8])*1e6,
                time=float.(data[:,5]),
                frequency=float.(data[:,4])*1e9,
                bandwidth=fill(4e6, size(data,1)),
                baselines=baselines
            )
    mjd =  Int(modified_julian(UTCEpoch(Int(data[1,2]), 1,1,0)).Δt+float(data[1,3]))
    return EHTObservation(data=edata, mjd=mjd, ra=180.0, dec=0.0, source=Symbol(data[1,1]))
end
