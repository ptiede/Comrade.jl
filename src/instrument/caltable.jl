export caltable

"""
    $(TYPEDEF)

A Tabes of calibration quantities. The columns of the table are the telescope sites codes.
The rows are the calibration quantities at a specific time and frequency. This user should not
use this struct directly. Instead that should call [`caltable`](@ref).
"""
struct CalTable{T,F,G<:AbstractVecOrMat}
    names::Vector{Symbol}
    lookup::Dict{Symbol,Int}
    times::T
    freqs::F
    gmat::G
end

function caltable end

Tables.istable(::Type{<:CalTable}) = true
Tables.columnaccess(::Type{<:CalTable}) = true

"""
    sites(g::CalTable)

Return the sites in the calibration table
"""
sites(g::CalTable) = getfield(g, :names)
times(g::CalTable) = getfield(g, :times)
frequencies(g::CalTable) = getfield(g, :freqs)
lookup(g::CalTable) = getfield(g, :lookup)
gmat(g::CalTable) = getfield(g, :gmat)
function Tables.schema(g::CalTable{T,F,G}) where {T,F,G}
    nms = [:Ti, :Fr]
    append!(nms, sites(g))
    types = Type[eltype(T), eltype(F)]
    append!(types, fill(eltype(G), size(gmat(g),2)))
    return Tables.Schema(nms, types)
end


Tables.columns(g::CalTable) = Tables.table(hcat(times(g), frequencies(g), gmat(g)); header=Tables.columnnames(g))
function Tables.getcolumn(g::CalTable, ::Type{T}, col::Int, nm::Symbol) where {T}
    (col == 1 || nm == :Ti) && return times(g)
    (col == 2 || nm == :Fr) && return frequencies(g)
    gmat(g)[:, col-2]
end

function Tables.getcolumn(g::CalTable, nm::Symbol)
    nm == :Ti && return times(g)
    nm == :Fr && return frequencies(g)
    return gmat(g)[:, lookup(g)[nm]]
end

function Tables.getcolumn(g::CalTable, i::Int)
    i==1 && return times(g)
    i==2 && return frequencies(g)
    return gmat(g)[:, i-2]
end

function viewcolumn(gt::CalTable, nm::Symbol)
    nm == :Ti && return times(gt)
    nm == :Fr && return frequencies(gt)
    return @view gmat(gt)[:, lookup(gt)[nm]]
end

Tables.columnnames(g::CalTable) = [:Ti, :Fr, sites(g)...]

Tables.rowaccess(::Type{<:CalTable}) = true
Tables.rows(g::CalTable) = g

struct CalTableRow{T,G} <: Tables.AbstractRow
    row::Int
    source::CalTable{T,G}
end

Tables.columnnames(g::CalTableRow) = [:Ti, :Fr, sites(getfield(g, :source))...]

function Base.propertynames(g::CalTable)
    return Tables.columnnames(g)
end

function Base.getproperty(g::CalTable, nm::Symbol)
    return viewcolumn(g, nm)
end

function Base.getindex(gt::CalTable, i::Int, nm::Symbol)
    nm == :Ti && return times(gt)[i]
    nm == :Fr && return frequencies(gt)[i]
    return Tables.getcolumn(gt, nm)[i]
end

function Base.getindex(gt::CalTable, ::typeof(Base.:!) , nm::Symbol)
    getproperty(gt, nm)
end

function Base.getindex(gt::CalTable, I::AbstractUnitRange, nm::Symbol)
    getproperty(gt, nm)[I]
end
function Base.getindex(gt::CalTable, I::AbstractVector{Int}, nm::Symbol)
    getproperty(gt, nm)[I]
end


function Base.view(gt::CalTable, I::AbstractUnitRange, nm::Symbol)
    @view getproperty(gt, nm)[I]
end

function Base.view(gt::CalTable, I::AbstractVector{Int}, nm::Symbol)
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
    (col == 1 || nm == :Ti) && return times(getfield(g, :source))[getfield(g, :row)]
    (col == 2 || nm == :Fr) && return frequencies(getfield(g, :source))[getfield(g, :row)]
    gmat(getfield(g, :source))[getfield(g, :row), col-2]
end

function Tables.getcolumn(g::CalTableRow, i::Int)
    (i==1) && return times(getfield(g, :source))[getfield(g, :row)]
    (i==2) && return frequencies(getfield(g, :source))[getfield(g, :row)]
    gmat(getfield(g, :source))[getfield(g, :row), i-2]
end

function Tables.getcolumn(g::CalTableRow, nm::Symbol)
    src = getfield(g, :source)
    nm == :Ti && return times(src)[getfield(g, :row)]
    nm == :Fr && return frequencies(src)[getfield(g, :row)]
    return gmat(src)[getfield(g, :row), lookup(src)[nm]]
end

@recipe function f(gt::CalTable; sites=Comrade.sites(gt), datagains=false)

    @argcheck prod(sites .∈ Ref(Comrade.sites(gt))) "passed site isn't in array\n"*
                                                "sites:     $(sites)\n"*
                                                "telescope: $(sites(gt))"
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
    t = getproperty.(gt[:Ti], :t0)
    xlims --> (t[begin], t[end] + 0.01*abs(t[end]))
    for (i,s) in enumerate(sites)
        T = nonmissingtype(eltype(gt[s]))
        for (j,f) in enumerate(unique(gt[:Fr]))
            @series begin
                seriestype := :scatter
                subplot := i
                title := String(s)
                if i == length(sites)
                    xguide --> "Time (UTC)"
                end
                label = "$(round(f.central/1e9, digits=1)) GHz"
                if i == length(sites)
                    label --> label
                else
                    label := nothing
                end
                ind = Base.:!.(ismissing.(gt[s]))
                find = findall(==(f), gt[:Fr][ind])
                #x := gt[:Ti][ind]
                x = t[ind][find]

                if !datagains
                    y = T.(gt[s][ind][find])
                else
                    y = T.(inv.((gt[s][ind][find])))
                end
                x, y
            end
        end
    end
end


using PrettyTables

function Base.show(io::IO, ct::CalTable, )
    pretty_table(io, Tables.columns(ct);
                     header=Tables.columnnames(ct),
                     vlines=[1,2],
                    formatters = _ctab_formatter
                )
end

function _ctab_formatter(v, i, j)
    j == 1 && return "$(round(v.t0, digits=2)) hr"
    j == 2 && return "$(round(v.central, digits=2)/1e9) GHz"
    return round(v, digits=3)
end


"""
    caltable(s::SiteArray)

Creates a calibration table from a site array
"""
function caltable(sarr::SiteArray)
    sites = sort(unique(Comrade.sites(sarr)))
    tf = collect(Iterators.product(unique(times(sarr))|>sort, unique(frequencies(sarr))|>sort))
    time = vec(first.(tf))
    freq = vec(last.(tf))
    gmat = Matrix{Union{eltype(sarr), Missing}}(missing, length(time), length(sites))
    gmat .= missing
    lookup = Dict(sites[i]=>i for i in eachindex(sites))
    for (j, s) in enumerate(sites)
        cterms = site(sarr, s)
        for (i, (t,f)) in enumerate(tf)
            ti = times(cterms)
            fi = frequencies(cterms)
            ind = findfirst(i->((ti[i]==t)&&fi[i]==f), eachindex(ti, fi))
            if !isnothing(ind)
                gmat[i, j] = cterms[ind]
            end
        end
    end
    return CalTable(sites, lookup, time, freq, gmat)
end
