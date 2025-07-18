export datatable

"""
    $(TYPEDEF)

An abstract VLBI table that is used to store data for a VLBI observation.
To implement your own table you just need to specify the  `VLBISkyModels.rebuild` function.
"""
abstract type AbstractVLBITable{F} end

"""
    datatable(obs::AbstractVLBITable)

Construct a table from the observation `obs`. The table is usually a StructArray of fields
"""
datatable(obs::AbstractVLBITable) = getfield(obs, :datatable)


function Base.getindex(config::F, i::AbstractVector) where {F <: AbstractVLBITable}
    newconf = datatable(config)[i]
    return rebuild(config, newconf)
end

function Base.view(config::F, i::AbstractVector) where {F <: AbstractVLBITable}
    newconf = @view(datatable(config)[i])
    return rebuild(config, newconf)
end


# Implement the tables interface
Tables.istable(::Type{<:AbstractVLBITable}) = true
Tables.columnaccess(::Type{<:AbstractVLBITable}) = true
Tables.columns(t::AbstractVLBITable) = datatable(t)

Tables.getcolumn(t::AbstractVLBITable, ::Type{T}, col::Int, nm::Symbol) where {T} = Tables.getcolumn(t, nm)
Tables.getcolumn(t::AbstractVLBITable, nm::Symbol) = getproperty(datatable(t), nm)
Tables.getcolumn(t::AbstractVLBITable, i::Int) = Tables.getcolumn(t, Tables.columnnames(t)[i])
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


abstract type AbstractMount end

"""
    Mount(parallactic, elevation, offset=0)

Defines the telescope mount type. The parallactic argument controls the evolution
of the parallactic angle, the elevation controls the elevation of the mount, 
and the offset is an optional offset of the feeds in radians.

For specific telescopes you can use the convenience constructors:
- `MountCassegrain(offset=0)`: for Cassegrain mounts
- `MountNasmythR(offset=0)`: for right Nasmyth mounts
- `MountNasymthL(offset=0)`: for left Nasmyth mounts
"""
struct Mount{A, B, C} <: AbstractMount
    parallactic::A
    elevation::B
    offset::C
    function Mount(parallactic::A, elevation::B, offset::C = 0.0) where {A, B, C}
        return new{A, B, C}(parallactic, elevation, offset)
    end
end

parallactic_mount(m::Mount) = getfield(m, :parallactic)
elevation_mount(m::Mount) = getfield(m, :elevation)
offset_mount(m::Mount) = getfield(m, :offset)

MountCassegrain(offset = 0) = Mount(1, 0, offset)
MountNasmythR(offset = 0) = Mount(1, 1, offset)
MountNasymthL(offset = 0) = Mount(1, -1, offset)


Base.@kwdef struct Antenna{N, T, M <: AbstractMount, F <: PolBasis}
    name::N
    position::SVector{3, T}
    mount::M
    feed::F
end

