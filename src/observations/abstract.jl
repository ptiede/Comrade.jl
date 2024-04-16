export datatable
abstract type AbstractVLBITable{F} end
datatable(obs::AbstractVLBITable) = getfield(obs, :datatable)

# Implement the tables interface
Tables.istable(::Type{<:AbstractVLBITable}) = true
Tables.columnaccess(::Type{<:AbstractVLBITable}) = true
Tables.columns(t::AbstractVLBITable) = getfield(t, :data)

Tables.getcolumn(t::AbstractVLBITable, ::Type{T}, col::Int, nm::Symbol) where {T} = getdata(t, nm)
Tables.getcolumn(t::AbstractVLBITable, nm::Symbol) = getproperty(datatable(t), nm)
Tables.getcolumn(t::AbstractVLBITable, i::Int) = Tables.getcolumn(t, Tables.columnames(t)[i])
Tables.columnnames(t::AbstractVLBITable) = propertynames(datatable(t))

Base.getindex(data::AbstractVLBITable, s::Symbol) = Tables.getcolumn(data, s)
Base.length(data::AbstractVLBITable) = length(datatable(data))
Base.lastindex(data::AbstractVLBITable) = lastindex(data.data)
Base.firstindex(data::AbstractVLBITable) = firstindex(data.data)
Base.getindex(data::AbstractVLBITable, i::Int) = build_datum(data, i)
