abstract type AbstractObservationTable{F<:AbstractVisibilityDatum} <: AbstractVLBITable{F} end
measurement(t::AbstractObservationTable) = getfield(t, :measurement)
noise(t::AbstractObservationTable) = getfield(t, :noise)
baseline(t::AbstractObservationTable) = datatable(arrayconfig(t))
arrayconfig(c::AbstractObservationTable, p::Symbol) = getindex(arrayconfig(c), p)
Base.length(obs::AbstractObservationTable) = length(measurement(obs))
Base.firstindex(obs::AbstractObservationTable) = firstindex(measurement(obs))
Base.lastindex(obs::AbstractObservationTable) = lastindex(measurement(obs))
datumtype(::AbstractObservationTable{T}) where {T} = T

function domain(obs::AbstractObservationTable)
    return domain(arrayconfig(obs))
end


function datatable(obs::AbstractObservationTable{F}) where {F}
    StructArray((build_datum(obs, i) for i in 1:length(obs)), unwrap=(T->(T<:Tuple || T<:AbstractBaselineDatum)))
end

"""
    arrayconfig(obs::AbstractObservationTable)

Returns the array configuration for a given observation.
"""
arrayconfig(obs::AbstractObservationTable) = getfield(obs, :config)

"""
    build_datum(data::AbstractObservationTable, i::Int)

Build the datum `F` for a given observation table `data`. This is
an internal method that users shouldn't have to deal with directly
unless they are implementing a new `AbstractObservationTable`.
"""
function build_datum(data::AbstractObservationTable{F}, i::Int) where {F}
    arr = arrayconfig(data)
    m   = measurement(data)
    e   = noise(data)
    return build_datum(F, m[i], e[i], arr[i])
end

"""
    beamsize(obs::AbstractObservationTable)

Calculate the approximate beam size of the observation `obs` as the inverse of the longest baseline
distance.
"""
beamsize(obs::AbstractObservationTable) = beamsize(arrayconfig(obs))

"""
    sites(d::AbstractObservationTable)

Get all the sites in a observation. The result is a vector of symbols.
"""
function sites(d::AbstractObservationTable)
    return sites(arrayconfig(d))
end



"""
    $(TYPEDEF)

The main data product type in `Comrade` this stores the `data` which can be a StructArray
of any `AbstractInterferometryDatum` type.

# Fields
$FIELDS
"""
struct EHTObservationTable{T<:AbstractVisibilityDatum,S<:AbstractArray, E<:AbstractArray, A<:AbstractArrayConfiguration} <: AbstractObservationTable{T}
    """
    Obervation measurement
    """
    measurement::S
    """
    Observation thermal noise
    """
    noise::E
    """
    Array config holds ancillary information about array
    """
    config::A
    function EHTObservationTable{T}(meas, err, config) where {T}
        return new{T, typeof(meas), typeof(err), typeof(config)}(meas, err, config)
    end
    function EHTObservationTable{A, B, C, D}(meas, err, config) where {A, B, C, D}
        return new{A, B, C, D}(meas, err, config)
    end
end

function VLBISkyModels.rebuild(data::EHTObservationTable{T}, newtable) where {T}
    m = newtable.measurement
    s = newtable.noise
    b = newtable.baseline
    newconf = rebuild(arrayconfig(data), b)
    return EHTObservationTable{T}(m, s, newconf)
end

function Base.show(io::IO, d::EHTObservationTable{F}) where {F}
    config = arrayconfig(d)
    sF = split("$F", ",")[1]
    sF = sF*"}"
    println(io, "EHTObservation{$sF}")
    println(io, "  source:      ", config.source)
    println(io, "  mjd:         ", config.mjd)
    println(io, "  bandwidth:   ", config.bandwidth)
    println(io, "  sites:       ", sites(d))
    print(io,   "  nsamples:    ", length(config))
end

function Base.getindex(obs::EHTObservationTable{F}, i::AbstractVector) where {F<:ClosureProducts}
    conf = arrayconfig(obs)[i]
    m = measurement(obs)[i]
    s = noise(obs)[i, i]
    return EHTObservationTable{F}(m, s, conf)
end

function Base.view(obs::EHTObservationTable{F}, i::AbstractVector) where {F<:ClosureProducts}
    conf = @view arrayconfig(obs)[i]
    m = @view measurement(obs)[i]
    s = @view noise(obs)[i, i]
    return EHTObservationTable{F}(m, s, conf)
end
