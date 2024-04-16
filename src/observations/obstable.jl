abstract type AbstractObservationTable{F<:AbstractVisibilityDatum} <: AbstractVLBITable{F} end
measurement(t::AbstractObservationTable) = getfield(t, :measurement)
error(t::AbstractObservationTable) = getfield(t, :error)
baseline(t::AbstractObservationTable) = datatable(arrayconfig(t))

function Base.getindex(data::F, i::AbstractVector) where {F<:AbstractObservationTable}
    config = arrayconfig(data)
    newconftable = datatable(config)[i]
    newconfig = @set config.datatable = newconftable
    return F(measurement(data)[i], error(data)[i], newconfig)
end


function datatable(obs::AbstractObservationTable{F}) where {F}
    return StructArray{F}((measurement=obs.measurement, error=obs.error, baseline=datatable(obs.config)))
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
    e   = error(data)
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
    Observation thermal error
    """
    error::E
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


function Base.show(io::IO, d::EHTObservationTable{F}) where {F}
    config = arrayconfig(d)
    sF = split("$F", ",")[1]
    sF = sF*"}"
    println(io, "EHTObservation{$sF}")
    println(io, "  source:      ", config.source)
    println(io, "  mjd:         ", config.mjd)
    println(io, "  frequencies: ", unique(config[:F]))
    println(io, "  bandwidth:   ", config.bandwidth)
    println(io, "  sites:       ", sites(config))
    print(io, "  nsamples:    ", length(config))
end
