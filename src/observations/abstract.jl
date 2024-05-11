export datatable
abstract type AbstractVLBITable{F} end
datatable(obs::AbstractVLBITable) = getfield(obs, :datatable)
function Base.getindex(config::F, i::AbstractVector) where {F<:AbstractVLBITable}
    newconf = datatable(config)[i]
    return rebuild(config, newconf)
end

function Base.view(config::F, i::AbstractVector) where {F<:AbstractVLBITable}
    newconf = @view(datatable(config)[i])
    return rebuild(config, newconf)
end


# Implement the tables interface
Tables.istable(::Type{<:AbstractVLBITable}) = true
Tables.columnaccess(::Type{<:AbstractVLBITable}) = true
Tables.columns(t::AbstractVLBITable) = getfield(t, :data)

Tables.getcolumn(t::AbstractVLBITable, ::Type{T}, col::Int, nm::Symbol) where {T} = getdata(t, nm)
Tables.getcolumn(t::AbstractVLBITable, nm::Symbol) = getproperty(datatable(t), nm)
Tables.getcolumn(t::AbstractVLBITable, i::Int) = Tables.getcolumn(t, Tables.columnames(t)[i])
Tables.columnnames(t::AbstractVLBITable) = propertynames(datatable(t))

function Base.getindex(data::AbstractVLBITable, s::Symbol)
    s == :measurement && return measurement(data)
    s == :noise       && return noise(data)
    return Tables.getcolumn(data, s)
end
Base.length(data::AbstractVLBITable) = length(datatable(data))
Base.lastindex(data::AbstractVLBITable) = lastindex(datatable(data))
Base.firstindex(data::AbstractVLBITable) = firstindex(datatable(data))
Base.getindex(data::AbstractVLBITable, i::Int) = build_datum(data, i)
function VLBISkyModels.rebuild(config::AbstractVLBITable, table)
    throw(MethodError(rebuild, config, table))
end
