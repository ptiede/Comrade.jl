export caltable

"""
    $(TYPEDEF)

A Tabes of calibration quantities. The columns of the table are the telescope station codes.
The rows are the calibration quantities at a specific time stamp. This user should not
use this struct directly. Instead that should call [`caltable`](@ref).
"""
struct CalTable{T,G<:AbstractVecOrMat}
    names::Vector{Symbol}
    lookup::Dict{Symbol,Int}
    times::T
    gmat::G
end

"""
    caltable(args...)

Creates a calibration table from a set of arguments. The specific arguments
depend on what calibration you are applying.
"""
function caltable end

Tables.istable(::Type{<:CalTable}) = true
Tables.columnaccess(::Type{<:CalTable}) = true

"""
    stations(g::CalTable)

Return the stations in the calibration table
"""
stations(g::CalTable) = getfield(g, :names)
scantimes(g::CalTable) = getfield(g, :times)
lookup(g::CalTable) = getfield(g, :lookup)
gmat(g::CalTable) = getfield(g, :gmat)
function Tables.schema(g::CalTable{T,G}) where {T,G}
    nms = [:time]
    append!(nms, stations(g))
    types = Type[eltype(T)]
    append!(types, fill(eltype(G), size(gmat(g),2)))
    return Tables.Schema(nms, types)
end


Tables.columns(g::CalTable) = Tables.table([scantimes(g) gmat(g)]; header=Tables.columnnames(g))
function Tables.getcolumn(g::CalTable, ::Type{T}, col::Int, nm::Symbol) where {T}
    (i == 1 || nm == :time) && return scantimes(g)
    gmat(g)[:, col-1]
end

function Tables.getcolumn(g::CalTable, nm::Symbol)
    nm == :time && return scantimes(g)
    return gmat(g)[:, lookup(g)[nm]]
end

function Tables.getcolumn(g::CalTable, i::Int)
    i==1 && return scantimes(g)
    return gmat(g)[:, i-1]
end

function viewcolumn(gt::CalTable, nm::Symbol)
    nm == :time && return scantimes(gt)
    return @view gmat(gt)[:, lookup(gt)[nm]]
end

Tables.columnnames(g::CalTable) = [:time, stations(g)...]

Tables.rowaccess(::Type{<:CalTable}) = true
Tables.rows(g::CalTable) = g

struct CalTableRow{T,G} <: Tables.AbstractRow
    row::Int
    source::CalTable{T,G}
end

function Base.propertynames(g::CalTable)
    return Tables.columnnames(g)
end

function Base.getproperty(g::CalTable, nm::Symbol)
    return viewcolumn(g, nm)
end

function Base.getindex(gt::CalTable, i::Int, nm::Symbol)
    nm == :time && return gt.times[i]
    return Tables.getcolumn(gt, nm)[i]
end

function Base.getindex(gt::CalTable, ::typeof(Base.:!) , nm::Symbol)
    getproperty(gt, nm)
end

function Base.getindex(gt::CalTable, I::AbstractUnitRange, nm::Symbol)
    getproperty(gt, nm)[I]
end

function Base.view(gt::CalTable, I::AbstractUnitRange, nm::Symbol)
    @view getproperty(gt, nm)[I]
end

function Base.getindex(gt::CalTable, ::Colon, nm::Symbol)
    Tables.getcolumn(gt, nm)
end

function Base.getindex(gt::CalTable, i::Int, ::Colon)
    Tables.getrow(gt, i)
end

Base.getindex(gt::CalTable, n::Symbol) = getproperty(gt, n)

Base.eltype(::CalTable{T,G}) where {T,G} = CalTableRow{T,G}
Base.length(g::CalTable) = size(gmat(g),1)
Base.iterate(m::CalTable, st=1) = st > length(m) ? nothing : (CalTableRow(st, m), st + 1)

function Tables.getrow(g::CalTable, i::Int)
    return CalTableRow(i, g)
end

function Tables.getrow(g::CalTable, I::AbstractUnitRange)
    [Tables.getrow(g, i) for i in I]
end


function Tables.getcolumn(g::CalTableRow, ::Type, col::Int, nm::Symbol)
    (col == 1 || nm == :time) && return scantimes(getfield(g, :source))[getfield(g, :row)]
    gmat(getfield(g, :source))[getfield(g, :row), col-1]
end

function Tables.getcolumn(g::CalTableRow, i::Int)
    (i==1) && return scantimes(getfield(g, :source))[getfield(g, :row)]
    gmat(getfield(g, :source))[getfield(g, :row), i-1]
end

function Tables.getcolumn(g::CalTableRow, nm::Symbol)
    src = getfield(g, :source)
    nm == :time && return scantimes(src)[getfield(g, :row)]
    return gmat(src)[getfield(g, :row), lookup(src)[nm]]
end

@recipe function f(gt::CalTable; sites=stations(gt), datagains=false)

    @argcheck prod(sites .∈ Ref(stations(gt))) "passed site isn't in array\n"*
                                                "sites:     $(sites)\n"*
                                                "telescope: $(stations(gt))"
    #if !datagains
    #    plot_title --> "Model Gain Amp."
    #else
    #    plot_title --> "Data Gain Amp."
    #end
    layout --> (length(sites), 1)


    size --> (350, 150*length(sites))
    #lims = extrema(filter(!ismissing, gmat(gt))).*1.1
    #if !datagains
    #    ylims --> lims
    #else
    #    ylims --> inv.(lims)[end:-1:begin]
    #end
    xlims --> (gt[:time][begin], gt[:time][end] + 0.01*abs(gt[:time][end]))
    for (i,s) in enumerate(sites)
        @series begin
            seriestype := :scatter
            subplot := i
            label --> :none

            if i == length(sites)
                xguide --> "Time (UTC)"
            end

            T = nonmissingtype(eltype(gt[s]))
            ind = Base.:!.(ismissing.(gt[s]))
            #x := gt[:time][ind]
            if !datagains
                yy = gt[s][ind]
            else
                yy = inv.(gt[s])[ind]
            end

            title --> string(s)
            gt[:time][ind], T.(yy)
        end
    end
end

Tables.columnnames(g::CalTableRow) = [:time, stations(getfield(g, :source))...]

using PrettyTables

function Base.show(io::IO, ct::CalTable, )
    pretty_table(io, Tables.columns(ct);
                     header=Tables.columnnames(ct),
                     vlines=[1],
                     formatters = (v,i,j)->round(v, digits=3)
                )
end


"""
    caltable(obs::EHTObservation, gains::AbstractVector)

Create a calibration table for the observations `obs` with `gains`. This returns
a [`CalTable`](@ref) object that satisfies the
`Tables.jl` interface. This table is very similar to the `DataFrames` interface.

# Example

```julia
ct = caltable(obs, gains)

# Access a particular station (here ALMA)
ct[:AA]
ct.AA

# Access a the first row
ct[1, :]
```
"""
function caltable(obs::EHTObservation, gains::AbstractVector, seg = ScanSeg(), segmented=false)
    gcache = jonescache(obs, seg, segmented)
    return caltable(gcache, gains)
end

"""
    caltable(g::JonesCache, jterms::AbstractVector)

Convert the JonesCache `g` and recovered Jones/corruption elements `jterms` into a `CalTable`
which satisfies the `Tables.jl` interface.

# Example

```julia
ct = caltable(gcache, gains)

# Access a particular station (here ALMA)
ct[:AA]
ct.AA

# Access a the first row
ct[1, :]
```
"""
function caltable(g::JonesCache, gains::AbstractVector)
    @argcheck length(g.times) == length(gains)

    stations = sort(unique(g.stations))
    times = unique(g.times)
    gmat = Matrix{Union{eltype(gains), Missing}}(missing, length(times), length(stations))

    alltimes = g.times
    allstations = g.stations
    # Create the lookup dict
    lookup = Dict(stations[i]=>i for i in eachindex(stations))
    for i in eachindex(gains)
        t = findfirst(t->(t==alltimes[i]), times)
        c = lookup[allstations[i]]
        gmat[t,c] = gains[i]
    end

    return CalTable(stations, lookup, times, gmat)
end

function caltable(g::SegmentedJonesCache, gains::AbstractVector)
    @argcheck length(g.times) == length(gains)

    stations = sort(unique(g.stations))
    times = unique(g.times)
    gmat = Matrix{Union{eltype(gains), Missing}}(undef, length(times), length(stations))
    gmat .= 0.0
    alltimes = g.times
    allstations = g.stations
    # Create the lookup dict
    lookup = Dict(stations[i]=>i for i in eachindex(stations))
    for i in eachindex(gains)
        t = findfirst(t->(t==alltimes[i]), times)
        c = lookup[allstations[i]]
        gmat[t,c] = gains[i]
    end


    cumsum!(gmat, gmat; dims=1)
    replace!(x->x==0 ? missing : x, @view gmat[:,begin:end])

    return CalTable(stations, lookup, times, gmat)
end
