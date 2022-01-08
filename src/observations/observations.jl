abstract type AbstractInterferometryDatum{T} end

abstract type AbstractVisibilityDatum{T} <: AbstractInterferometryDatum{T} end
abstract type AbstractLinearPolDatum{S<:AbstractVisibilityDatum, T} <: AbstractInterferometryDatum{T} end
abstract type AbstractCrossPolDatum{S,T} <: AbstractInterferometryDatum{T} end

abstract type ClosureProducts{T} <: AbstractInterferometryDatum{T} end

abstract type Observation{T} end

using DelimitedFiles
using AstroTime: modified_julian

export uvpositions, stations, getdata, nsamples, arrayconfig,
        nsamples, getuv, baselines



"""
    $(TYPEDEF)
This defined the abstract type for an array configuration. Namely, baseline
times, SEFD's, bandwidth, observation frequencies, etc.
"""
abstract type ArrayConfiguration end


struct EHTArrayConfiguration{S,F,T<:AbstractArray} <: ArrayConfiguration
    stations::S
    frequency::F
    bandwidth::F
    uvsamples::T
end

function getuv(ac::ArrayConfiguration)
    return ac.uvsamples.u, ac.uvsamples.v
end

function uvtimefreq(ac::EHTArrayConfiguration)
    u,v = getuv(ac)
    t = ac.uvsamples.time
    ν = ac.frequency
    return u, v, t, fill(ν, length(u))
end

struct ArrayBaselineDatum{T}
    time::T
    u::T
    v::T
    baseline::Tuple{Symbol, Symbol}
    error_real::T
    error_imag::T
end

const ArrayQuadrangleDatum = NTuple{4, ArrayBaselineDatum{T}} where {T}
const ArrayTriangleDatum = NTuple{3, ArrayBaselineDatum{T}} where {T}


uvpositions(D::AbstractVisibilityDatum) = D.u, D.v

Base.@kwdef struct EHTObservation{F,T<:AbstractInterferometryDatum{F},S<:StructArray{T}, N} <: Observation{F}
    data::S
    mjd::N
    ra::F
    dec::F
    bandwidth::F
    frequency::F
    source::Symbol
    timetype::Symbol = :UTC
end

Base.getindex(data::EHTObservation, i::Int) = data.data[i]
Base.length(data::EHTObservation) = nsamples(data)

function stations(d::EHTObservation{T,A}) where {T,A<:AbstractInterferometryDatum}
    bl = getdata(d, :baseline)
    s1 = first.(bl)
    s2 = last.(bl)
    return unique([s1..., s2...])
end


function Base.show(io::IO, d::EHTObservation{F,D}) where {F,D}
    println(io, "EHTObservation{$F,$D, ...}")
    println(io, "  source: ", d.source)
    println(io, "  mjd: ", d.mjd)
    println(io, "  frequency: ", d.frequency)
    println(io, "  bandwidth: ", d.bandwidth)
    println(io, "  stations: ", stations(d))
    println(io, "  nsamples: ", nsamples(d))
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
    baseline::NTuple{2,Symbol}
end

function amplitude(D::EHTVisibilityDatum)
    amp = hypot(D.visr,D.visi)
    return EHTVisibilityAmplitudeDatum(amp, D.error,
                                       D.u, D.v,
                                       D.time,
                                       D.frequency,
                                       D.bandwidth,
                                       D.baseline
                                    )
end

Base.@kwdef struct EHTVisibilityAmplitudeDatum{T<:Number} <: AbstractVisibilityDatum{T}
    amp::T
    error::T
    u::T
    v::T
    time::T
    baseline::NTuple{2,Symbol}
end


function checktriangle(D1::EHTVisibilityDatum,
                       D2::EHTVisibilityDatum,
                       D3::EHTVisibilityDatum)
    b1 = D1.baseline
    b2 = D2.baseline
    b3 = D3.baseline
    l = length(unique([b1..., b2..., b3...]))
    @assert l == 3 "For a valid closure phase you need 3 unique stations not $l"
    @assert (D1.time == D2.time == D3.time) "For a valid closure phase the times need to match"

end

@inline function visibility(D::EHTVisibilityDatum{T}) where {T}
    return Complex{T}(D.visr, D.visi)
end

@inline function amplitude(D::EHTVisibilityAmplitudeDatum{T}) where {T}
    return D.amp
end

@inline function bispectrum(D1::EHTVisibilityDatum, D2::EHTVisibilityDatum, D3::EHTVisibilityDatum)
    visibility(D1)*visibility(D2)*visibility(D3)
end

"""
    $(TYPEDEF)

A Datum for a single closure phase. Note in the future this may get replaced
with the fully covariant formalism from Blackburn et al. (2020).
"""
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
    triangle::NTuple{3,Symbol}
end

function stations(d::EHTObservation{T,A}) where {T,A<:EHTClosurePhaseDatum}
    bl = getdata(d, :triangle)
    return unique(vcat(collect.(bl)...))
end




"""
    $(SIGNATURES)
Computes the closure phase of the three visibility datums.

# Notes
We currently use the high SNR Gaussian error approximation for the closure phase.
In the future we may use the moment matching from Monte Carlo sampling.
"""
function closure_phase(D1::EHTVisibilityDatum,
                      D2::EHTVisibilityDatum,
                      D3::EHTVisibilityDatum)

    checktriangle(D1,D2,D3)

    amp1 = amplitude(D1).amp
    amp2 = amplitude(D2).amp
    amp3 = amplitude(D3).amp
    u1,v1 = uvpositions(D1)
    u2,v2 = uvpositions(D2)
    u3,v3 = uvpositions(D3)

    bis = bispectrum(D1, D2, D3)
    s12 = unique([D1.baseline..., D2.baseline...])
    s123 = unique([s12..., D3.baseline...])
    #Use the Gaussian approximation TODO hook this into Measurements.jl?
    err = sqrt((D1.error/amp1)^2 + (D2.error/amp2)^2 + (D3.error/amp3)^2)
    return EHTClosurePhaseDatum(angle(bis), err,
                                u1, v1, u2, v2, u3, v3,
                                time, s123)
end

function baselines(CP::EHTClosurePhaseDatum)
    tri = CP.triangle
    return ((tri[1],tri[2]), (tri[2], tri[3]), (tri[3], tri[1]))
end

uvpositions(datum::EHTClosurePhaseDatum) = (datum.u1, datum.v1, datum.u2, datum.v2, datum.u3, datum.v3)


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
    quadrangle::NTuple{4,Symbol}
end

function baselines(CP::EHTLogClosureAmplitudeDatum)
    quad = CP.quadrangle
    return ((quad[1],quad[2]), (quad[3], quad[4]), (quad[1], quad[3]), (quad[2], quad[4]))
end

function stations(d::EHTObservation{T,A}) where {T,A<:EHTLogClosureAmplitudeDatum}
    bl = getdata(d, :quadrangle)
    return unique(vcat(collect.(bl)...))
end

uvpositions(datum::EHTLogClosureAmplitudeDatum) = (datum.u1, datum.v1, datum.u2, datum.v2, datum.u3, datum.v3, datum.u4, datum.v4)


# @doc raw"""
#     $(SIGNATURES)
# Compute the log-closure amplitude of 4 visibility datums.
# This uses the same ordering as eht-imaging i.e.

# ```math
# C_A = \log \frac{|V_{1}| |V_{2}|}{|V_{3}| |V_{4}|}
# ```

# # Notes

# Currently we use the high SNR Gaussian approximation to the thermal error.
# In the future we may use the moment matching approach using Monte Carlo sampling.
# """
# function logclosure_amplitude(D1::EHTVisibilityDatum,
#                               D2::EHTVisibilityDatum,
#                               D3::EHTVisibilityDatum,
#                               D4::EHTVisibilityDatum)
#     check_quadrangle(D1,D2,D3,D4)

#     amp1 = amplitude(D1)
#     amp2 = amplitude(D2)
#     amp3 = amplitude(D3)
#     amp4 = amplitude(D4)

#     u1,v1 = uvpositions(D1)
#     u2,v2 = uvpositions(D2)
#     u3,v3 = uvpositions(D3)
#     u4,v4 = uvpositions(D4)

#     #Construct the quadrangle
#     s12 = unique([D1.baseline..., D2.baseline...])
#     s123 = unique([s12..., D3.baseline...])
#     s1234 = unique([s123..., D4.baseline...])

#     # Now do the error propogation
#     err = hypot(D1.err/map1, D2.err/amp2, D3.err/amp3, D4.err/amp4)

#     #construct the lcamp
#     lcamp = log(amp1*amp2/(amp3*amp4))

#     return EHTLogClosureAmplitudeDatum(lcamp, err,
#                                        u1, v1,
#                                        u2, v2,
#                                        u3, v3,
#                                        u4, v4,
#                                        time, s1234)
# end


function arrayconfig(vis::EHTObservation{F,A}) where {F,A<:Union{EHTVisibilityDatum, EHTVisibilityAmplitudeDatum}}
    u = getdata(vis, :u)
    v = getdata(vis, :v)
    st = stations(vis)
    times = getdata(vis, :time)
    bandwidth = vis.bandwidth
    frequency = vis.frequency
    error = getdata(vis, :error)
    baseline = getdata(vis, :baseline)
    uvsamples = StructArray{ArrayBaselineDatum}(time=times,
                                        u=u,
                                        v=v,
                                        baseline=baseline,
                                        error_real=error,
                                        error_imag=error
                                    )
    return EHTArrayConfiguration(st, frequency, bandwidth, uvsamples)
end

include(joinpath(@__DIR__, "io.jl"))
include(joinpath(@__DIR__, "ehtim.jl"))
