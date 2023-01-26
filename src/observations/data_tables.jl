abstract type DataTables{T<:AbstractInterferometryDatum} end

scans(d::DataTables)           = scans(arrayconfig(d))
telescope_array(d::DataTables) = telescope_array(arrayconfig(d))
bandwidth(d::DataTables)       = bandwidth(arrayconfig(d))
measurement(d::DataTables)     = d.measurement
error(d::DataTables)           = d.error
arrayconfig(d::DataTables)     = d.config
obsdata(d::DataTables)         = d.obs
ra(d::DataTables)              = d.ra
dec(d::DataTables)             = d.dec
source(d::DataTables)          = d.source
stations(d::EHTDataTable) = stations(arrayconfig(d))


"""
    data(obs::DataTables)
    data(obs::DataTables, s::Symbol)

Pass-through function that gets the array of `s` from the DataTables.
If no `s` is supplied it returns a StructArray of datums that correspond to the
specific data product stored in the table.

For example say you want the times of all measurement then
```julia
data(obs, :time)
```
"""
function data(obs::DataTables{F})
    build.(F, arrayconfig(obs), measurement(obs), error(obs))
end

function data(obs::DataTables, s::Symbol)
    s === :measurement && return measurement(obs)
    s === :error       && return error(obs)
    return data(arrayconfig(obs), s)
end

getuv(table::DataTables) = getuv(arrayconfig(table))

Base.length(data::DataTables) = length(measurement(data))
Base.eltype(::DataTables{T}) where T = T


# Implement the tables interface
Tables.istable(::Type{<:DataTables}) = true
Tables.columnaccess(::Type{<:DataTables}) = true
Tables.columns(t::DataTables) = Tables.table(data(t))

Tables.getcolumn(t::DataTables, ::Type{T}, col::Int, nm::Symbol) where {T} = data(t, nm)
Tables.getcolumn(t::DataTables, nm::Symbol) = data(t, nm)
Tables.getcolumn(t::DataTables, i::Int) = Tables.getcolumn(t, Tables.columnames(t)[i])
Tables.columnnames(t::DataTables) = propertynames(data(t))

Tables.rowaccess(::Type{<:DataTables}) = true
Tables.rows(t::DataTable) = t

function Tables.getrow(table::DataTables{F}, i::Int) where F
    c = arrayconfig(table)[i]
    m = measurement(table)[i]
    e = error(table)[i]
    return build(F, c, m, e)
end

Tables.getrow(table::DataTables{F}, i::AbstractUnitRange) = Tables.getrow.(table, i)
Base.getindex(t::DataTables, i::Int, ::Colon) = Tables.getrow(t, i)


function Base.show(io::IO, ct::DataTables)
    pretty_table(io, Tables.columns(ct);
                     header=Tables.columnnames(ct),
                     vlines=[1],
                     formatters = (v,i,j)->round(v, digits=3),
                     title = "  source:    $(source(ct))\n"*
                             "  mjd:       $(mjd(ct))\n"*
                             "  bandwidth: $(bandwidth(ct))\n"
                )
end
