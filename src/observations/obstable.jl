export domain, datatable, arrayconfig, sites, beamsize, EHTObservationTable, snr, uvdist,
    measwnoise

"""
    $(TYPEDEF)

The abstract obervation table. This contains actual data plus the array configuration.
"""
abstract type AbstractObservationTable{F <: AbstractVisibilityDatum} <: AbstractVLBITable{F} end
measurement(t::AbstractObservationTable) = getfield(t, :measurement)
noise(t::AbstractObservationTable) = getfield(t, :noise)
baseline(t::AbstractObservationTable) = datatable(arrayconfig(t))


uvdist(d) = hypot(d.baseline.U, d.baseline.V)

function uvdist(d::EHTClosurePhaseDatum)
    u = map(x -> x.U, d.baseline)
    v = map(x -> x.V, d.baseline)
    a = hypot(u[1] - u[2], v[1] - v[2])
    b = hypot(u[2] - u[3], v[2] - v[3])
    c = hypot(u[3] - u[1], v[3] - v[1])
    return sqrt(heron(a, b, c))
end

function heron(a, b, c)
    s = 0.5 * (a + b + c)
    return sqrt(s * (s - a) * (s - b) * (s - c))
end

function uvdist(d::EHTLogClosureAmplitudeDatum)
    u = map(x -> x.U, d.baseline)
    v = map(x -> x.V, d.baseline)
    a = hypot(u[1] - u[2], v[1] - v[2])
    b = hypot(u[2] - u[3], v[2] - v[3])
    c = hypot(u[3] - u[4], v[3] - v[4])
    d = hypot(u[4] - u[1], v[4] - v[1])
    h = hypot(u[1] - u[3], v[1] - v[3])
    return sqrt(heron(a, b, h) + heron(c, d, h))
end


arrayconfig(c::AbstractObservationTable, p::Symbol) = getindex(arrayconfig(c), p)
Base.length(obs::AbstractObservationTable) = length(measurement(obs))
Base.firstindex(obs::AbstractObservationTable) = firstindex(measurement(obs))
Base.lastindex(obs::AbstractObservationTable) = lastindex(measurement(obs))

# Returns the data type of the observation table
datumtype(::AbstractObservationTable{T}) where {T} = T


"""
    domain(obs::AbstractObservationTable; executor=Serial(), header=ComradeBase.NoHeader()

Returns the u, v, time, frequency domain of the observation.
"""
function domain(obs::AbstractObservationTable; executor = Serial(), header = ComradeBase.NoHeader())
    return domain(arrayconfig(obs); executor, header)
end

"""
    datatable(obs::AbstractObservationTable)

Returns a tabular representation of the data. Note that for closures this ignores the covariance
between quantities, which is otherwise included in the full `EHTObservationTable`.
"""
function datatable(obs::AbstractObservationTable{F}) where {F}
    return StructArray((build_datum(obs, i) for i in 1:length(obs)), unwrap = (T -> (T <: Tuple || T <: AbstractBaselineDatum || T <: SArray || T <: NamedTuple)))
end

"""
    arrayconfig(obs::AbstractObservationTable)
    arrayconfig(obs::AbstractObservationTable, p::Symbol)

Returns the array configuration for a given observation. If `p` is provided then only the
property `p` is returned.
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
    m = measurement(data)
    e = noise(data)
    return build_datum(F, m[i], e[i], arr[i])
end

function build_datum(data::AbstractObservationTable{F}, i::Int) where {F <: ClosureProducts}
    arr = arrayconfig(data)
    m = measurement(data)
    e = noise(data)
    return build_datum(F, m[i], sqrt(e[i, i]), arr[i])
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
of any `AbstractInterferometryDatum` type. Note that the underlying structure is not part
of the public API. Users should typically construct tables from the [`extract_table`](@ref)
function.

# Fields
$FIELDS
"""
struct EHTObservationTable{T <: AbstractVisibilityDatum, S <: AbstractArray, E <: AbstractArray, A <: AbstractArrayConfiguration} <: AbstractObservationTable{T}
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
    return EHTObservationTable{eltype(newtable)}(m, s, newconf)
end

function Base.show(io::IO, d::EHTObservationTable{F}) where {F}
    config = arrayconfig(d)
    sF = split("$F", ",")[1]
    sF = sF * "}"
    println(io, "EHTObservationTable{$sF}")
    println(io, "  source:      ", config.source)
    println(io, "  mjd:         ", config.mjd)
    println(io, "  bandwidth:   ", config.bandwidth)
    println(io, "  sites:       ", sites(d))
    return print(io, "  nsamples:    ", length(config))
end

function Base.getindex(obs::EHTObservationTable{F}, i::AbstractVector) where {F <: ClosureProducts}
    conf = arrayconfig(obs)[i]
    m = measurement(obs)[i]
    s = noise(obs)[i, i]
    return EHTObservationTable{F}(m, s, conf)
end

function Base.view(obs::EHTObservationTable{F}, i::AbstractVector) where {F <: ClosureProducts}
    conf = @view arrayconfig(obs)[i]
    m = @view measurement(obs)[i]
    s = @view noise(obs)[i, i]
    return EHTObservationTable{F}(m, s, conf)
end

"""
    reset_mounts!(d::EHTObservation, overrides)

Apply per-site overrides to the array configuration of an `EHTObservation` from a dict of the form `Dict{Symbol, NamedTuple}` 
where the keys are site names and the values are named tuples with fields `SEFD1`, `SEFD2`, `fr_parallactic`, `fr_elevation`, and `fr_offset`.
You can construct the overrides from ehtim-style array.txt files via [`Comrade.load_array_txt`](@ref) or programmatically.
"""
function reset_mounts!(d::EHTObservationTable, overrides)
    tarr = arrayconfig(d).tarr
    for i in eachindex(tarr.sites)
        site = tarr.sites[i]
        if haskey(overrides, site)
            o = overrides[site]
            tarr.SEFD1[i] = o.SEFD1
            tarr.SEFD2[i] = o.SEFD2
            tarr.fr_parallactic[i] = o.fr_parallactic
            tarr.fr_elevation[i] = o.fr_elevation
            tarr.fr_offset[i] = o.fr_offset
        end
    end
    return d
end
