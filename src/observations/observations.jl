abstract type AbstractInterferometryDatum{T} end

abstract type AbstractVisibilityDatum{T} <: AbstractInterferometryDatum{T} end
abstract type AbstractLinearPolDatum{S<:AbstractVisibilityDatum, T} <: AbstractInterferometryDatum{T} end
abstract type AbstractCrossPolDatum{S,T} <: AbstractInterferometryDatum{T} end

abstract type ClosureProducts{T} <: AbstractInterferometryDatum{T} end

abstract type Observation{T} end

using DelimitedFiles
using AstroTime: UTCEpoch, modified_julian


"""
    $(TYPEDEF)
This defined the abstract type for an array configuration. Namely, baselines
times, SEFD's, bandwidth, observation frequencies, etc.
"""
abstract type ArrayConfiguration end


struct StaticArray{S} <: ArrayConfiguration
    sourcemeta::S
    arraydata::D
end


getuv(D::AbstractVisibilityDatum) = D.u, D.v

Base.@kwdef struct EHTObservation{F,T<:AbstractInterferometryDatum{F},S<:StructArray{T}, N} <: Observation{F}
    data::S
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

vamp(d::EHTVisibilityDatum) = hypot(d.visr, d.visi)
phase(d::EHTVisibilityDatum) = angle(d.visr + im*d.visi)

function getamps(D::EHTVisibilityDatum)
    amp = hypot(D.visr,D.visi)
    return EHTVisibilityAmplitudeDatum(amp, D.error,
                                       D.u, D.v,
                                       D.time,
                                       D.frequency,
                                       D.bandwidth,
                                       D.baselines
                                    )
end

Base.@kwdef struct EHTVisibilityAmplitudeDatum{T<:Number} <: AbstractVisibilityDatum{T}
    amp::T
    error::T
    u::T
    v::T
    time::T
    frequency::T
    bandwidth::T
    baselines::NTuple{2,Symbol}
end

#=
function closurephase(D1::EHTVisibilityDatum,
                      D2::EHTVisibilityDatum,
                      D3::EHTVisibilityDatum)

    checktriangle(D1,D2,D3)

    u1,v2 = getuv(D1)
    amp1 = vamp(D1)
    amp2 = vamp(D2)
    amp3 = vamp(D3)
    u2,v2 = getuv(D1)
    u3,v3 = getuv(D1)

    bis = bispectrum(D1,D2,D3)
    err = sqrt( (D1.error/amp1)^2 + (D2.error/amp2)^2 + (D3.error/amp3)^2)
    return EHTClosurePhaseDatum(angle(bis), err, u1,v1,u2,v2,u3,v3,
end
=#

Base.@kwdef struct EHTLogClosureAmplitudeDatum{T<:Number} <: ClosureProducts{T}
    amp::T
    error::T
    u1::T
    v1::T
    u2::T
    v2::T
    u3::T
    v3::T
    u4::T
    v4::T
    time::T
    frequency::T
    bandwidth::T
    baselines::NTuple{4,Symbol}
end



Base.@kwdef struct EHTClosurePhaseDatum{T<:Number} <: ClosureProducts{T}
    phase::T
    error::T
    u1::T
    v1::T
    u2::T
    v2::T
    u3::T
    v3::T
    time::T
    frequency::T
    bandwidth::T
    baselines::NTuple{3,Symbol}
end
uvtriangle(datum::EHTClosurePhaseDatum) = (datum.u1, datum.v1, datum.u2, datum.v2, datum.u3, datum.v3)

"""
    $(SIGNATURES)
Load a ThemisPy style ascii EHT observation file.
"""
function load_tpy(file)
    data = readdlm(file, skipstart=1)
    bs1 = Symbol.(getindex.(data[:,5], Ref(1:2)))
    bs2 = Symbol.(getindex.(data[:,5], Ref(3:4)))
    baselines = tuple.(bs1, bs2)
    edata = StructArray{EHTVisibilityDatum{Float64}}(
                visr=float.(data[:,8]),
                visi=float.(data[:,10]),
                error=float.(data[:,9]),
                u=float.(data[:,6])*1e6,
                v=float.(data[:,7])*1e6,
                time=float.(data[:,4]),
                frequency=fill(227e9, size(data,1)),
                bandwidth=fill(4e6, size(data,1)),
                baselines=baselines
            )
    mjd =  Int(modified_julian(UTCEpoch(Int(data[1,2]), 1,1,0)).Δt+float(data[1,3]))
    return EHTObservation(data=edata, mjd=mjd, ra=180.0, dec=0.0, source=Symbol(data[1,1]))
end
