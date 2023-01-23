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
       getuv, baselines, scantable



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
struct EHTArrayConfiguration{F,T,D<:AbstractArray} <: ArrayConfiguration
    """
    Observing bandwith (Hz)
    """
    bandwidth::F
    """
    Telescope array file
    """
    tarr::T
    """
    A struct array of `ArrayBaselineDatum` holding time, freq, u, v, baselines.
    """
    data::D
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
function getuv(ac::ArrayConfiguration)
    return (U=ac.data.U, V=ac.data.V)
end



function getuv(ac::ClosureConfig)
    return getuv(ac.ac)
end

"""
    $(SIGNATURES)

Get the u, v, time, freq of the array as a tuple.
"""
function getuvtimefreq(ac::EHTArrayConfiguration)
    u,v = getuv(ac)
    t = ac.data.T
    ν = ac.data.F
    return (U=u, V=v, T=t, F=ν)
end

function getuvtimefreq(ac::ClosureConfig)
    return getuvtimefreq(ac.ac.config)
end


"""
    $(TYPEDEF)

A single datum of an `ArrayConfiguration`
"""
struct ArrayBaselineDatum{T,E,V}
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
    """
    elevation of baselines
    """
    elevation::Tuple{V,V}
    """
    parallactic angle of baslines
    """
    parallactic::Tuple{V,V}
    function ArrayBaselineDatum(u, v, time, freq, baseline, error, elevation, parallactic)
        tt, ft, ut, vt = promote(time, freq, u, v)
        T = typeof(tt)
        V = typeof(elevation[1])
        E = typeof(error)
        return new{T,E,V}(ut, vt, tt, ft, baseline, error, elevation, parallactic)
    end
end



const ArrayQuadrangleDatum = NTuple{4, ArrayBaselineDatum{T}} where {T}
const ArrayTriangleDatum = NTuple{3, ArrayBaselineDatum{T}} where {T}

"""
    uvpositions(datum::AbstractVisibilityDatum)

Get the uvp positions of an inferometric datum.
"""
uvpositions(D::AbstractVisibilityDatum) = D.U, D.V

"""
    $(TYPEDEF)

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

function getuv(ac::EHTObservation)
    return (U=ac.data.U, V=ac.data.V)
end

function getuv(ac::EHTObservation{T,A}) where {T,A<:ClosureProducts}
    return (U=ac.config.ac.data.U, V=ac.config.ac.data.V)
end


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
    return sort(unique(vcat(s1, s2)))
end


function Base.show(io::IO, d::EHTObservation{F,D}) where {F,D}
    println(io, "EHTObservation{$F,$D, ...}")
    println(io, "  source: ", d.source)
    println(io, "  mjd: ", d.mjd)
    println(io, "  frequency: ", first(d.data.F))
    println(io, "  bandwidth: ", d.bandwidth)
    println(io, "  stations: ", stations(d))
    println(io, "  nsamples: ", length(d))
end






"""
    $(TYPEDEF)

A struct holding the information for a single measured visibility.

# $(FIELDS)

"""
Base.@kwdef struct EHTVisibilityDatum{S<:Number} <: AbstractVisibilityDatum{S}
    """
    real component of the visibility (Jy)
    """
    measurement::Complex{S}
    """
    error of the visibility (Jy)
    """
    error::S
    """
    u position of the data point in λ
    """
    U::S
    """
    v position of the data point in λ
    """
    V::S
    """
    time of the data point in (Hr)
    """
    T::S
    """
    frequency of the data point (Hz)
    """
    F::S
    """
    station baseline codes
    """
    baseline::NTuple{2,Symbol}
end

"""
    visibility(d::EHTVisibilityDatum)

Return the complex visibility of the visibility datum
"""
@inline function visibility(D::EHTVisibilityDatum{T}) where {T}
        return D.measurement
end


"""
    amplitude(d::EHTVisibilityDatum)

Get the amplitude of a visibility datum
"""
function amplitude(D::EHTVisibilityDatum)
    amp = hypot(D.visr,D.visi)
    return EHTVisibilityAmplitudeDatum(amp, D.error,
                                       D.u, D.v,
                                       D.time,
                                       D.frequency,
                                       D.baseline
                                    )
end

"""
    $(TYPEDEF)

A struct holding the information for a single measured visibility amplitude.

# FIELDS
$(FIELDS)

"""
Base.@kwdef struct EHTVisibilityAmplitudeDatum{S<:Number} <: AbstractVisibilityDatum{S}
    """
    amplitude (Jy)
    """
    measurement::S
    """
    error of the visibility amplitude (Jy)
    """
    error::S
    """
    u position of the data point in λ
    """
    U::S
    """
    v position of the data point in λ
    """
    V::S
    """
    time of the data point in (Hr)
    """
    T::S
    """
    frequency of the data point (Hz)
    """
    F::S
    """
    station baseline codes
    """
    baseline::NTuple{2,Symbol}
end

# internal method that checks whether the triangle is closes
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

"""
    amplitude(d::EHTVisibilityAmplitudeDatum)

Get the amplitude of a amplitude datum
"""
@inline function amplitude(D::EHTVisibilityAmplitudeDatum{T}) where {T}
    return D.measurement
end


"""
    bispectrum(d1::T, d2::T, d3::T) where {T<:EHTVisibilityDatum}

Finds the bispectrum of three visibilities. We will assume these form closed triangles,
i.e. the phase of the bispectrum is a closure phase.
"""
@inline function bispectrum(D1::EHTVisibilityDatum, D2::EHTVisibilityDatum, D3::EHTVisibilityDatum)
    checktriangle(D1, D2, D3)
    visibility(D1)*visibility(D2)*visibility(D3)
end



"""
    $(TYPEDEF)

A Datum for a single closure phase.

# Fields
$(FIELDS)

"""
Base.@kwdef struct EHTClosurePhaseDatum{S<:Number} <: ClosureProducts{S}
    """
    closure phase (rad)
    """
    measurement::S
    """
    error of the closure phase assuming the high-snr limit
    """
    error::S
    """
    u (λ) of first station
    """
    U1::S
    """
    v (λ) of first station
    """
    V1::S
    """
    u (λ) of second station
    """
    U2::S
    """
    v (λ) of second station
    """
    V2::S
    """
    u (λ) of third station
    """
    U3::S
    """
    v (λ) of third station
    """
    V3::S
    """
    Measured time of closure phase in hours
    """
    T::S
    """
    Measured frequency of closure phase in Hz
    """
    F::S
    """
    station baselines used
    """
    triangle::NTuple{3,Symbol}
end


"""
    $(TYPEDEF)

A Datum for a single coherency matrix

# Fields
$(FIELDS)

"""
Base.@kwdef struct EHTCoherencyDatum{S, B1, B2, M<:SMatrix{2,2,Complex{S}}, E<:SMatrix{2,2,S}} <: Comrade.AbstractInterferometryDatum{S}
    """
    coherency matrix, with entries in Jy
    """
    measurement::M
    """
    visibility uncertainty matrix, with entries in Jy
    """
    error::E
    """
    x-direction baseline length, in λ
    """
    U::S
    """
    y-direction baseline length, in λ
    """
    V::S
    """
    Timestamp, in hours
    """
    T::S
    """
    Frequency, in Hz
    """
    F::S
    """
    station baseline codes
    """
    baseline::NTuple{2,Symbol}
    """
    polarization basis for each station
    """
    polbasis::Tuple{B1, B2}
end


function stations(d::EHTObservation{T,A}) where {T,A<:EHTClosurePhaseDatum}
    bl = getdata(d, :triangle)
    return sort(unique(vcat(collect.(bl)...)))
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
Base.@kwdef struct EHTLogClosureAmplitudeDatum{S<:Number} <: ClosureProducts{S}
    """
    log-closure amplitude
    """
    measurement::S
    """
    log-closure amplitude error in the high-snr limit
    """
    error::S
    """
    u (λ) of first station
    """
    U1::S
    """
    v (λ) of first station
    """
    V1::S
    """
    u (λ) of second station
    """
    U2::S
    """
    v (λ) of second station
    """
    V2::S
    """
    u (λ) of third station
    """
    U3::S
    """
    v (λ) of third station
    """
    V3::S
    """
    u (λ) of fourth station
    """
    U4::S
    """
    v (λ) of fourth station
    """
    V4::S
    """
    Measured time of closure phase in hours
    """
    T::S
    """
    Measured frequency of closure phase in Hz
    """
    F::S
    """
    station codes for the quadrangle
    """
    quadrangle::NTuple{4,Symbol}
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
    return sort(unique(vcat(collect.(bl)...)))
end

uvpositions(datum::EHTLogClosureAmplitudeDatum) = (datum.U1, datum.V1, datum.U2, datum.V2, datum.U3, datum.V3, datum.U4, datum.V4)


"""
    $(SIGNATURES)

Extract the array configuration from a EHT observation.
"""
function arrayconfig(vis::EHTObservation)
    vis.config
end


function _arrayconfig(data, angles, tarr, bandwidth)
    u = getproperty(data, :U)
    v = getproperty(data, :V)
    times = getproperty(data, :T)
    error = getproperty(data, :error)
    baseline = getproperty(data, :baseline)
    frequency = getproperty(data, :F)
    uvsamples = StructArray{ArrayBaselineDatum}(T=times,
                                        U=u,
                                        V=v,
                                        F = frequency,
                                        baseline=baseline,
                                        error=error,
                                        elevation = StructArray(angles[1]),
                                        parallactic  = StructArray(angles[2])
                                    )
    return EHTArrayConfiguration(bandwidth, tarr, uvsamples)
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
function baselines(scan::Scan{A,B,C}) where {A,B,C<:StructArray{<:AbstractInterferometryDatum}}
    bl = scan.scan.baseline
    # organize the closure phase stations
    ant1 = first.(bl)
    ant2 = last.(bl)
    return ant1, ant2
end


# Closures are special
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




function stations(s::Scan)
    ants = baselines(s)
    stat = unique(vcat(ants...))
    return sort(stat)
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
    times = obs[:T]
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
