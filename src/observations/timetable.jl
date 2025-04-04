export timetable

"""
    $(TYPEDEF)

Wraps EHTObservationTable in a table that separates the observation into scans.
This implements the table interface. You can access scans by directly
indexing into the table. This will create a view into the table not copying the data.

# Example
```julia-repl
julia> st = timetable(obs)
julia> st[begin] # grab first scan
julia> st[end]   # grab last scan
```
"""
struct TimeTable{O <: AbstractVLBITable, T, S}
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

Base.length(st::TimeTable) = length(st.times)
Base.firstindex(st::TimeTable) = firstindex(st.times)
Base.lastindex(st::TimeTable) = lastindex(st.times)
sites(st::TimeTable) = sites(st.obs)
Base.eachindex(st::TimeTable) = LinearIndices(firstindex(st):lastindex(st))
function Base.getindex(st::TimeTable, i::Int)
    istart = st.scanind[i]
    if i < length(st.scanind)
        iend = st.scanind[i + 1] - 1
    else
        iend = length(st.obs)
    end
    return Scan(st.times[i], (i, istart, iend), @view st.obs[istart:iend])
end


"""
    $(TYPEDEF)

Composite type that holds information for a single scan of the telescope.

# Fields
$(FIELDS)
"""
struct Scan{T, I, S}
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
    baseline(scan::Scan)

Return the baselines for each datum in a scan
"""
function baseline(scan::Scan{A, B, <:EHTObservationTable}) where {A, B}
    bl = arrayconfig(scan.scan)[:sites]
    # organize the closure phase sites
    ant1 = first.(bl)
    ant2 = last.(bl)
    return ant1, ant2
end


function sites(s::Scan)
    return sites(s.scan)
end

function Base.show(io::IO, s::Scan)
    println(io, "VLBI Scan")
    println(io, "\tscan index: ", s.index)
    println(io, "\tscan time:  ", s.time)
    return print(io, "\tsites: ", sites(s))
end


function Base.getindex(st::TimeTable, I)
    return [getindex(st, i) for i in I]
end

function Base.getindex(scan::Scan, s::Symbol)
    return getindex(scan.scan, s)
end

function Base.getindex(scan::Scan, i::Int)
    return scan.scan[i]
end

function Base.getindex(scan::Scan, i::AbstractVector{<:Union{Bool, Int}})
    return Scan(scan.time, scan.index, scan.scan[i])
end


Base.first(st::TimeTable) = st[1]
Base.last(st::TimeTable) = st[length(st)]

"""
    timetable(obs::AbstractArrayConfiguration)
Reorganizes the observation into a table of scans, where scan are defined by unique timestamps.
To access the data you can use scalar indexing

# Example

```julia
st = timetable(obs)
# Grab the first scan
scan1 = st[1]

# Acess the detections in the scan
scan1[1]

# grab e.g. the baselines
scan1[:baseline]
```

"""
function timetable(arr::AbstractVLBITable)
    times = domain(arr).Ti
    scantimes = unique(times)
    scanind = Int[]
    for t in scantimes
        ind = findfirst(==(t), times)
        append!(scanind, ind)
    end
    return TimeTable(arr, scantimes, scanind)
end
