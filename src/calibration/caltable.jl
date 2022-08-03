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

Tables.columns(g::CalTable) = [scantimes(g) gmat(g)]
Tables.getcolumn(g::CalTable, ::Type{T}, col::Int, nm::Symbol) where {T} = gmat(g)[:, col]
function Tables.getcolumn(g::CalTable, nm::Symbol)
    nm == :time && return scantimes(g)
    return gmat(g)[:, lookup(g)[nm]]
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

function Base.getindex(gt, I::AbstractUnitRange, nm::Symbol)
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
    gmat(getfield(g, :source))[getfield(m, :row), col]
end

function Tables.getcolumn(g::CalTableRow, i::Int)
    gmat(getfield(g, :source))[getfield(g, :row), i]
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
    for (i,s) in enumerate(sites)
        @series begin
            seriestype := :scatter
            subplot := i
            label := :none

            if i == length(sites)
                xguide --> "Time (UTC)"
            end

            T = nonmissingtype(eltype(gt[s]))
            ind = Base.:!.(ismissing.(gt[s]))
            x := gt[:time][ind]
            if !datagains
                yy = gt[s][ind]
            else
                yy = inv.(gt[s])[ind]
            end

            title --> string(s)
            T.(yy)
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
