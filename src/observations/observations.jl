"""
    $(TYPEDEF)
An abstract type for all VLBI interfermetry data types. See [EHTVisibilityDatum](@ref) for an example.
"""
abstract type AbstractInterferometryDatum{T} end

abstract type AbstractVisibilityDatum{T} <: AbstractInterferometryDatum{T} end
abstract type AbstractLinearPolDatum{S<:AbstractVisibilityDatum, T} <: AbstractInterferometryDatum{T} end
abstract type AbstractCrossPolDatum{S,T} <: AbstractInterferometryDatum{T} end

abstract type ClosureProducts{T} <: AbstractInterferometryDatum{T} end

abstract type Observation{T} end

using DelimitedFiles
using AstroTime: modified_julian

export uvpositions, stations, getdata, arrayconfig,
       getuv, baselines, rescaleuv!, scantable



"""
    $(TYPEDEF)

This defined the abstract type for an array configuration. Namely, baseline
times, SEFD's, bandwidth, observation frequencies, etc.
"""
abstract type ArrayConfiguration end

"""
    $(TYPEDEF)

Stores all the non-visibility data products for an EHT array. This is useful when evaluating
model visibilities.

# Fields
$(FIELDS)
"""
struct EHTArrayConfiguration{F,T<:AbstractArray} <: ArrayConfiguration
    """
    Observing frequency (Hz)
    """
    frequency::F
    """
    Observing bandwith (Hz)
    """
    bandwidth::F
    """
    A struct array of `ArrayBaselineDatum` holding time, freq, u, v, baselines.
    """
    data::T
end

"""
    $(TYPEDEF)
Array config file for closure quantities. This stores the design matrix `designmat`
that transforms from visibilties to closure products.
"""
struct ClosureConfig{A,D} <: ArrayConfiguration
    ac::A
    designmat::D
function ClosureConfig(ac, dmat)
    A = typeof(ac)
    sdmat = blockdiag(sparse.(dmat)...)
    D = typeof(sdmat)
    return new{A,D}(ac, sdmat)
end
end


"""
    $(SIGNATURES)
Get the u, v positions of the array.
"""
function getuv(ac::ArrayConfiguration)
    return ac.data.u, ac.data.v
end

function getuv(ac::ClosureConfig)
    return getuv(ac.ac)
end

"""
    $(SIGNATURES)
Get the u, v, time, freq of the array as a tuple.
"""
function uvtimefreq(ac::EHTArrayConfiguration)
    u,v = getuv(ac)
    t = ac.data.time
    ν = ac.frequency
    return u, v, t, fill(ν, length(u))
end

"""
    $(TYPEDEF)

A single datum of an `ArrayConfiguration`
"""
struct ArrayBaselineDatum{T}
    """
    time of the data point in (Hr)
    """
    time::T
    """
    frequency of the data point (Hz)
    """
    freq::T
    """
    u position of the data point in λ
    """
    u::T
    """
    v position of the data point in λ
    """
    v::T
    """
    Station codes of the baseline (u,v)
    """
    baseline::Tuple{Symbol, Symbol}
    function ArrayBaselineDatum(time, freq, u, v, baseline)
        tt, ft, ut, vt = promote(time, freq, u, v)
        T = typeof(tt)
        return new{T}(tt, ft, ut, vt, baseline)
    end
end



const ArrayQuadrangleDatum = NTuple{4, ArrayBaselineDatum{T}} where {T}
const ArrayTriangleDatum = NTuple{3, ArrayBaselineDatum{T}} where {T}

"""
    uvpositions(datum)
Get the uvp positions of an inferometric datum.
"""
uvpositions(D::AbstractVisibilityDatum) = D.u, D.v

"""
    $(SIGNATURES)
The main data product type in `Comrade` this stores the `data` which can be a StructArray
of any `AbstractInterferometryDatum` type.

# Fields
$FIELDS
"""
Base.@kwdef struct EHTObservation{F,T<:AbstractInterferometryDatum{F},S<:StructArray{T}, A, N} <: Observation{F}
    """
    StructArray of data productts
    """
    data::S
    """
    Array config holds ancillary information about array
    """
    config::A
    """
    modified julia date of the observation
    """
    mjd::N
    """
    RA of the observation in J2000 (deg)
    """
    ra::F
    """
    DEC of the observation in J2000 (deg)
    """
    dec::F
    """
    bandwidth of the observation (Hz)
    """
    bandwidth::F
    """
    frequency of the observation (Hz)
    """
    frequency::F
    """
    Common source name
    """
    source::Symbol
    """
    Time zone used.
    """
    timetype::Symbol = :UTC
end

getdata(obs::Observation, s::Symbol) = getproperty(getfield(obs, :data), s)


# Implement the tables interface
Tables.istable(::Type{<:EHTObservation}) = true
Tables.columnaccess(::Type{<:EHTObservation}) = true
Tables.columns(t::EHTObservation) = getfield(t, :data)

Tables.getcolumn(t::EHTObservation, ::Type{T}, col::Int, nm::Symbol) where {T} = getdata(t, nm)
Tables.getcolumn(t::EHTObservation, nm::Symbol) = getdata(t, nm)
Tables.getcolumn(t::EHTObservation, i::Int) = Tables.getcolumn(t, Tables.columnames(t)[i])
Tables.columnnames(t::EHTObservation) = propertynames(getfield(t, :data))

Base.getindex(data::EHTObservation, s::Symbol) = Tables.getcolumn(data, s)
Base.getindex(data::EHTObservation, i::Int) = data.data[i]
Base.getindex(data::EHTObservation, I...) = getindex(data.data, I...)
Base.length(data::EHTObservation) = length(data.data)

"""
    $(SIGNATURES)
Get all the stations in a observation. The result is a vector of symbols.
"""
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
    println(io, "  nsamples: ", length(d))
end

"""
    getdata(obs::EHTObservation, s::Symbol)
Pass-through function that gets the array of `s` from the EHTObservation. For example
say you want the times of all measurement then

```julia
getdata(obs, :time)
```
"""





"""
    $(SIGNATURES)
A struct holding the information for a single measured visibility.

# $(FIELDS)

"""
Base.@kwdef struct EHTVisibilityDatum{T<:Number} <: AbstractVisibilityDatum{T}
    """
    real component of the visibility (Jy)
    """
    visr::T
    """
    imaginary component of the visibility (Jy)
    """
    visi::T
    """
    error of the visibility (Jy)
    """
    error::T
    """
    x-direction baseline length in λ
    """
    u::T
    """
    y-direction baseline length in λ
    """
    v::T
    """
    Time of the observation in hours
    """
    time::T
    """
    station baseline codes
    """
    baseline::NTuple{2,Symbol}
end

"""
    $(SIGNATURES)
Return the complex visibility of the visibility datum
"""
@inline function visibility(D::EHTVisibilityDatum{T}) where {T}
        return Complex{T}(D.visr, D.visi)
end


"""
    $(SIGNATURES)
Get the amplitude of a visibility datum
"""
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

"""
    $(SIGNATURES)
A struct holding the information for a single measured visibility amplitude.

# $(FIELDS)

"""
Base.@kwdef struct EHTVisibilityAmplitudeDatum{T<:Number} <: AbstractVisibilityDatum{T}
    """
    amplitude (Jy)
    """
    amp::T
    """
    error of the visibility amplitude (Jy)
    """
    error::T
    """
    x-direction baseline length in λ
    """
    u::T
    """
    y-direction baseline length in λ
    """
    v::T
    """
    Time of the observation in hours
    """
    time::T
    """
    station baseline codes
    """
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


@inline function amplitude(D::EHTVisibilityAmplitudeDatum{T}) where {T}
    return D.amp
end

@inline function bispectrum(D1::EHTVisibilityDatum, D2::EHTVisibilityDatum, D3::EHTVisibilityDatum)
    visibility(D1)*visibility(D2)*visibility(D3)
end


"""
    rescaleuv!(data)
rescale the u-v lengths according to scale. This can be useful when you want the spatial scales
to be in 1/μas instead of 1/rad
"""
function rescaleuv!(data::EHTObservation{T,D}, scale) where {T,D<:Union{EHTVisibilityAmplitudeDatum,EHTVisibilityDatum}}
    data.data.u .= data.data.u.*scale
    data.data.v .= data.data.v.*scale
    return data
end



"""
    $(TYPEDEF)

A Datum for a single closure phase.

# $(FIELDS)

"""
Base.@kwdef struct EHTClosurePhaseDatum{T<:Number} <: ClosureProducts{T}
    """
    closure phase (rad)
    """
    phase::T
    """
    error of the closure phase assuming the high-snr limit
    """
    error::T
    """
    u (λ) of first station
    """
    u1::T
    """
    v (λ) of first station
    """
    v1::T
    """
    u (λ) of second station
    """
    u2::T
    """
    v (λ) of second station
    """
    v2::T
    """
    u (λ) of third station
    """
    u3::T
    """
    v (λ) of third station
    """
    v3::T
    """
    Measured time of closure phase in hours
    """
    time::T
    """
    station baselines used
    """
    triangle::NTuple{3,Symbol}
end

function rescaleuv!(data::EHTObservation{T,D}, scale) where {T,D<:EHTClosurePhaseDatum}
    data.data.u1 .= data.data.u1.*scale
    data.data.v1 .= data.data.v1.*scale
    data.data.u2 .= data.data.u2.*scale
    data.data.v2 .= data.data.v2.*scale
    data.data.u3 .= data.data.u3.*scale
    data.data.v3 .= data.data.v3.*scale
    return data
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

"""
    $(TYPEDEF)

A Datum for a single log closure amplitude.

# $(FIELDS)

"""
Base.@kwdef struct EHTLogClosureAmplitudeDatum{T<:Number} <: ClosureProducts{T}
    """
    log-closure amplitude
    """
    amp::T
    """
    log-closure amplitude error in the high-snr limit
    """
    error::T
    """
    u (λ) of first station
    """
    u1::T
    """
    v (λ) of first station
    """
    v1::T
    """
    u (λ) of second station
    """
    u2::T
    """
    v (λ) of second station
    """
    v2::T
    """
    u (λ) of third station
    """
    u3::T
    """
    v (λ) of third station
    """
    v3::T
    """
    u (λ) of fourth station
    """
    u4::T
    """
    v (λ) of fourth station
    """
    v4::T
    """
    Observation time of the quadrangle
    """
    time::T
    """
    station codes for the quadrangle
    """
    quadrangle::NTuple{4,Symbol}
end

function rescaleuv!(data::EHTObservation{T,D}, scale) where {T,D<:EHTLogClosureAmplitudeDatum}
    data.data.u1 .= data.data.u1.*scale
    data.data.v1 .= data.data.v1.*scale
    data.data.u2 .= data.data.u2.*scale
    data.data.v2 .= data.data.v2.*scale
    data.data.u3 .= data.data.u3.*scale
    data.data.v3 .= data.data.v3.*scale
    data.data.u4 .= data.data.u4.*scale
    data.data.v4 .= data.data.v4.*scale
    return data
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

"""
    $(SIGNATURES)
Extract the array configuration from a visibility EHT observation.
"""
function arrayconfig(vis::EHTObservation)
    vis.config
end

function _arrayconfig(data, bandwidth, frequency)
    u = getproperty(data, :u)
    v = getproperty(data, :v)
    times = getproperty(data, :time)
    error = getproperty(data, :error)
    baseline = getproperty(data, :baseline)
    uvsamples = StructArray{ArrayBaselineDatum}(time=times,
                                        u=u,
                                        v=v,
                                        freq = fill(frequency, length(u)),
                                        baseline=baseline,
                                        error_real=error,
                                        error_imag=error
                                    )
    return EHTArrayConfiguration(frequency, bandwidth, uvsamples)
end

const VisAmpDatum = Union{EHTVisibilityAmplitudeDatum, EHTVisibilityDatum}

"""
    `$(TYPEDEF)`
Wraps EHTObservation in a table that separates the observation into scans.
This implements the table interface and is used primarly for keeping track of
gains and closure information.
"""
struct ScanTable{O<:Union{Observation,ArrayConfiguration}, T, S}
    obs::O
    times::T
    scanind::S
end

Base.length(st::ScanTable) = length(st.times)
stations(st::ScanTable) = stations(st.obs)

struct Scan{T,I,S}
    time::T
    index::I
    scan::S
end

Base.length(s::Scan) = length(s.scan)

function baselines(scancp::Scan{A,B,C}) where {A,B,C<:StructArray{<:EHTClosurePhaseDatum}}
    tri = scancp.scan.triangle
    # organize the closure phase stations
    ant1 = getindex.(tri, 1)
    ant2 = getindex.(tri, 2)
    ant3 = getindex.(tri, 3)
    return ant1, ant2, ant3
end

function baselines(scancp::Scan{A,B,C}) where {A,B,C<:StructArray{<:EHTLogClosureAmplitudeDatum}}
    tri = scancp.scan.quadrangle
    # organize the closure phase stations
    ant1 = getindex.(tri, 1)
    ant2 = getindex.(tri, 2)
    ant3 = getindex.(tri, 3)
    ant4 = getindex.(tri, 4)
    return ant1, ant2, ant3, ant4
end

function baselines(scancp::Scan{A,B,C}) where {A,B,C<:StructArray{<:EHTVisibilityDatum}}
    bl = scancp.scan.baseline
    # organize the closure phase stations
    ant1 = first.(bl)
    ant2 = last.(bl)
    return ant1, ant2
end

function baselines(scancp::Scan{A,B,C}) where {A,B,C<:StructArray{<:EHTVisibilityAmplitudeDatum}}
    bl = scancp.scan.baseline
    # organize the closure phase stations
    ant1 = first.(bl)
    ant2 = last.(bl)
    return ant1, ant2
end


function stations(s::Scan)
    ants = baselines(s)
    stat = unique(vcat(ants...))
    return stat
end

function Base.show(io::IO, s::Scan)
    println(io, "VLBI Scan")
    println(io, "\tscan index: ", s.index)
    println(io, "\tscan time:  ", s.time)
    println(io, "\tstations: ", stations(s))
end

function Base.getindex(st::ScanTable, i::Int)
    istart = st.scanind[i]
    if i < length(st.scanind)
        iend = st.scanind[i+1]-1
    else
        iend = length(st.obs)
    end
    return Scan(st.times[i], (i, istart, iend), @view st.obs.data[istart:iend])
end

function Base.getindex(st::ScanTable, I)
    [getindex(st, i) for i in I]
end

function Base.getindex(scan::Scan, s::Symbol)
    getproperty(scan.scan, s)
end

function Base.getindex(scan::Scan, i::Int)
    scan.scan[i]
end

function Base.getindex(scan::Scan, i::AbstractVector{<:Union{Bool,Int}})
    Scan(scan.time, scan.index, scan.scan[i])
end


Base.first(st::ScanTable) = st[1]

"""
    `scantable(obs)`
Reorganizes the observation into a table of scans, where scan are defined by unique timestamps.
To access the data you can use scalar indexing

# Example

```julia
st = scantable(obs)
# Grab the first scan
scan1 = st[1]

# Acess the detections in the scan
scan1[1]

# grab e.g. the baselines
scan1[:baseline]
```

"""
function scantable(obs)
    times = obs[:time]
    scantimes = unique(times)
    scanind = Int[]
    for t in scantimes
        ind = findfirst(==(t), times)
        append!(scanind, ind)
    end
    return ScanTable(obs, scantimes, scanind)
end


include(joinpath(@__DIR__, "io.jl"))
include(joinpath(@__DIR__, "ehtim.jl"))
