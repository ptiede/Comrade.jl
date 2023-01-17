"""
    $(TYPEDEF)
An abstract type for all VLBI interfermetry data types. See [EHTVisibilityDatum](@ref) for an example.
"""
abstract type AbstractInterferometryDatum{T} end

abstract type AbstractVisibilityDatum{T} <: AbstractInterferometryDatum{T} end
abstract type ClosureProducts{T} <: AbstractInterferometryDatum{T} end

measurment(d::AbstractVisibilityDatum) = d.measurement
error(d::AbstractVisibilityDatum) = d.error
baseline(d::AbstractVisibilityDatum) = d.baseline
position(d::AbstractVisibilityDatum) = (U=d.U, V=d.V, T=d.T, F=d.F)
uvposition(d::AbstractVisibilityDatum) = (U=d.U, V=d.V)
feed(d::AbstractVisibilityDatum) = d.polbasis

abstract type Observation{T} end

using DelimitedFiles
using AstroTime: modified_julian

export uvpositions, stations, getdata, arrayconfig,
       getuv, baselines, rescaleuv!, scantable




struct EHTObservation
    telescope::T
    configuration::C
    measurements::M
    header::H
end


struct EHTArrayConfiguration
    U::TU
    V::TU
    T::TU
    F::TU
end

struct EHTHeader
    source::Symbol
    ra::F
    dec::F
    mjd::Int
end



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
    Observing bandwith (Hz)
    """
    bandwidth::F
    """
    A struct array of `ArrayBaselineDatum` holding time, freq, u, v, baselines.
    """
    data::T
end

function getindex(ac::EHTArrayConfiguration, v::Symbol)
    return getproperty(ac.data, v)
end

function getindex(ac::EHTArrayConfiguration, I::Union{Int, <:AbstractVector{<:Integer}})
    return getindex(ac, I)
end

"""
    $(TYPEDEF)
Array config file for closure quantities. This stores the design matrix `designmat`
that transforms from visibilties to closure products.

# Fields
$(FIELDS)
"""
struct ClosureConfig{A,D} <: ArrayConfiguration
    """Array configuration for visibilities"""
    ac::A
    """Closure design matrix"""
    designmat::D
function ClosureConfig(ac, dmat)
    A = typeof(ac)
    sdmat = blockdiag(sparse.(dmat)...)
    D = typeof(sdmat)
    return new{A,D}(ac, sdmat)
end
end


"""
    getuv

Get the u, v positions of the array.
"""
function uvposition(ac::ArrayConfiguration)
    return (U=ac.data.U, V=ac.data.V)
end



function uvposition(ac::ClosureConfig)
    return uvposition(ac.ac)
end

"""
    $(SIGNATURES)

Get the u, v, time, freq of the array as a tuple.
"""
function position(ac::EHTArrayConfiguration)
    (;U,V) = getuv(ac)
    t = ac.data.time
    ν = ac.frequency
    return (U=U, V=V, T=t, F=fill(ν, length(u)))
end

"""
    $(TYPEDEF)

A single datum of an `ArrayConfiguration`
"""
struct ArrayBaselineDatum{T, E} <: AbstractInterferometryDatum{T}
    """
    u position of the data point in λ
    """
    U::T
    """
    v position of the data point in λ
    """
    V::T
    """
    time of the data point in (Hr)
    """
    T::T
    """
    frequency of the data point (Hz)
    """
    F::T
    """
    Station codes of the baseline (u,v)
    """
    baseline::Tuple{Symbol, Symbol}
    """
    The thermal noise on the baseline
    """
    error::E
    function ArrayBaselineDatum(time, freq, u, v, baseline, error)
        tt, ft, ut, vt, errort = promote(time, freq, u, v, error)
        T = typeof(tt)
        return new{T}(tt, ft, ut, vt, baseline, errort)
    end
end


@non_differentiable getuv(ac::ArrayConfiguration)


const ArrayQuadrangleDatum = NTuple{4, ArrayBaselineDatum{T}} where {T}
const ArrayTriangleDatum = NTuple{3, ArrayBaselineDatum{T}} where {T}


"""
    $(TYPEDEF)

The main data product type in `Comrade` this stores the `data` which can be a StructArray
of any `AbstractInterferometryDatum` type.

# Fields
$FIELDS
"""
Base.@kwdef struct EHTObservation{F,T<:AbstractInterferometryDatum{F},S, A, N} <: Observation{F}
    """
    StructArray of data productts
    """
    measurement::S
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
    Common source name
    """
    source::Symbol
    """
    Time zone used.
    """
    timetype::Symbol = :UTC
end

"""
    getdata(obs::EHTObservation, s::Symbol)

Pass-through function that gets the array of `s` from the EHTObservation. For example
say you want the times of all measurement then

```julia
getdata(obs, :time)
```
"""
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
    stations(d::EHTObservation)

Get all the stations in a observation. The result is a vector of symbols.
"""
function stations(d::EHTObservation{T,A}) where {T,A<:AbstractInterferometryDatum}
    bl = getdata(d, :baseline)
    s1 = first.(bl)
    s2 = last.(bl)
    return unique(vcat(s1, s2))
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



abstract type MeasurementDatum{T} <: AbstractVisibilityDatum{T} end

error(d::MeasurementDatum) = error(d.baseline)
baseline(d::MeasurementDatum) = baseline(d.baseline)
position(d::MeasurementDatum) = position(d.baseline)
uvposition(d::MeasurementDatum) = uvposition(d.baseline)
feed(d::MeasurementDatum) = feed(d.baseline)


"""
    $(TYPEDEF)

A struct holding the information for a single measured visibility.

# $(FIELDS)

"""
Base.@kwdef struct VisibilityDatum{T<:Number, B<:ArrayBaselineDatum{T}} <: MeasurementDatum{T}
    """
    Complex visibility
    """
    measurement::Complex{T}
    """
    Baseline datum
    """
    config::B
end


"""
    $(TYPEDEF)

A struct holding the information for a single measured visibility amplitude.

# FIELDS
$(FIELDS)

"""
Base.@kwdef struct AmplitudeDatum{T<:Number} <: MeasurementDatum{T}
    """
    amplitude (Jy)
    """
    measurement::T
    """
    Baseline datum
    """
    config::B
end

# internal method that checks whether the triangle is closes
function checktriangle(D1::VisibilityDatum,
                       D2::VisibilityDatum,
                       D3::VisibilityDatum)
    b1 = baseline(D1)
    b2 = baseline(D2)
    b3 = baseline(D3)
    l = length(unique([b1..., b2..., b3...]))
    @assert l == 3 "For a valid closure phase you need 3 unique stations not $l"
    @assert (D1.time == D2.time == D3.time) "For a valid closure phase the times need to match"

end

"""
    amplitude(d::EHTVisibilityAmplitudeDatum)

Get the amplitude of a amplitude datum
"""
@inline function amplitude(D::VisibilityDatum{T}) where {T}
    return abs(measurement(D))
end


"""
    bispectrum(d1::T, d2::T, d3::T) where {T<:EHTVisibilityDatum}

Finds the bispectrum of three visibilities. We will assume these form closed triangles,
i.e. the phase of the bispectrum is a closure phase.
"""
@inline function bispectrum(D1::EHTVisibilityDatum, D2::EHTVisibilityDatum, D3::EHTVisibilityDatum)
    checktriangle(D1, D2, D3)
    measurement(D1)*measurement(D2)*measurement(D3)
end


"""
    rescaleuv!(data::EHTObservation)

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

# Fields
$(FIELDS)

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
    closure_phase(D1::EHTVisibilityDatum,
                  D2::EHTVisibilityDatum,
                  D3::EHTVisibilityDatum
                  )

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

"""
    baselines(CP::EHTClosurePhaseDatum)

Returns the baselines used for a single closure phase datum
"""
function baselines(CP::EHTClosurePhaseDatum)
    tri = CP.triangle
    return ((tri[1],tri[2]), (tri[2], tri[3]), (tri[3], tri[1]))
end

uvpositions(datum::EHTClosurePhaseDatum) = (datum.U1, datum.V1, datum.U2, datum.V2, datum.U3, datum.V3)

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

"""
    baselines(CP::EHTLogClosureAmplitudeDatum)

Returns the baselines used for a single closure phase datum
"""
function baselines(CP::EHTLogClosureAmplitudeDatum)
    quad = CP.quadrangle
    return ((quad[1],quad[2]), (quad[3], quad[4]), (quad[1], quad[3]), (quad[2], quad[4]))
end

function stations(d::EHTObservation{T,A}) where {T,A<:EHTLogClosureAmplitudeDatum}
    bl = getdata(d, :quadrangle)
    return unique(vcat(collect.(bl)...))
end

uvpositions(datum::EHTLogClosureAmplitudeDatum) = (datum.U1, datum.V1, datum.U2, datum.V2, datum.U3, datum.V3, datum.U4, datum.V4)


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

Extract the array configuration from a EHT observation.
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
    uvsamples = StructArray{ArrayBaselineDatum}(T=times,
                                        U=u,
                                        V=v,
                                        F = fill(frequency, length(u)),
                                        baseline=baseline,
                                        error=error
                                    )
    return EHTArrayConfiguration(frequency, bandwidth, uvsamples)
end

const VisAmpDatum = Union{EHTVisibilityAmplitudeDatum, EHTVisibilityDatum}

"""
    $(TYPEDEF)

Wraps EHTObservation in a table that separates the observation into scans.
This implements the table interface. You can access scans by directly
indexing into the table. This will create a view into the table not copying the data.

# Example
```julia-repl
julia> st = scantable(obs)
julia> st[begin] # grab first scan
julia> st[end]   # grab last scan
```
"""
struct ScanTable{O<:Union{Observation,ArrayConfiguration}, T, S}
    """
    Parent array information
    """
    obs::O
    """
    Scan times
    """
    times::T
    """
    Scan indices
    """
    scanind::S
end

Base.length(st::ScanTable) = length(st.times)
Base.firstindex(st::ScanTable) = firstindex(st.times)
Base.lastindex(st::ScanTable) = lastindex(st.times)
stations(st::ScanTable) = stations(st.obs)

"""
    $(TYPEDEF)

Composite type that holds information for a single scan of the telescope.

# Fields
$(FIELDS)
"""
struct Scan{T,I,S}
    """
    Scan time
    """
    time::T
    """
    Scan indices which are (scan index, data start index, data end index)
    """
    index::I
    """
    Scan data usually a StructArray of a <:AbstractVisibilityDatum
    """
    scan::S
end

Base.length(s::Scan) = length(s.scan)

"""
    baselines(scan::Scan)

Return the baselines for each datum in a scan
"""
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
Base.last(st::ScanTable) = st[length(st)]

"""
    scantable(obs::EHTObservation)
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
function scantable(obs::EHTObservation)
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