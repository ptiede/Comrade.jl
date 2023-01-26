"""
    $(TYPEDEF)

The main data product type in `Comrade` this stores the `data` which can be a StructArray
of any `AbstractInterferometryDatum` type.

# Fields
$FIELDS
"""
Base.@kwdef struct EHTDataTable{F,T<:AbstractInterferometryDatum{F},S,M,A,N,O<:Observation} <: DataTable{T}
    """
    Vector of measurements of the telescope
    """
    measurement::S
    """
    An array the holds the array of measurement errors
    """
    error::M
    """
    Array config holds ancillary information about array
    """
    config::A
    """
    RA of the observation in J2000 (deg)
    """
    ra::F
    """
    DEC of the observation in J2000 (deg)
    """
    dec::F
    """
    Common source name
    """
    source::Symbol
    """
    Base observation object
    """
    obs::O
end



"""
    stations(d::EHTDataTable)

Get all the stations in a observation. The result is a vector of symbols.
"""
function stations(d::EHTDataTable{T,A}) where {T,A<:AbstractInterferometryDatum}
    bl = data(d, :baseline)
    s1 = first.(bl)
    s2 = last.(bl)
    return sort(unique(vcat(s1, s2)))
end


function stations(d::EHTDataTable{T,A}) where {T,A<:EHTClosurePhaseDatum}
    bl = getdata(d, :triangle)
    return sort(unique(vcat(collect.(bl)...)))
end

function stations(d::EHTDataTable{T,A}) where {T,A<:EHTLogClosureAmplitudeDatum}
    bl = getdata(d, :quadrangle)
    return sort(unique(vcat(collect.(bl)...)))
end

"""
    scantable(obs::EHTDataTable)
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
function scantable(obs::EHTDataTable, dt=:scan)
    dt === :scan && return scantable_full(obs)


    tstart = first(times):dt:(last(times[end-1]))
    tstop  = times[begin+1]:dt:last(times[end]+dt)
    times = obs[:T]
    scanind = Int[]
    for t in scantimes
        ind = findfirst(t->tstart[i] <= t < tstop, times)
        append!(scanind, ind)
    end
    return ScanTable(obs, scantimes, scanind)
end


function scantable_full(obs::EHTDataTable)
    times = obs[:T]
    scantimes = scans(obs)
    tstart = scantimes.start
    tstop = scantimes.tstop
    scanind = Int[]
    for t in scantimes
        ind = findfirst(t->tstart[i] <= t <= tstop[i] , times)
        append!(scanind, ind)
    end
    return ScanTable(obs, scantimes, scanind)
end
