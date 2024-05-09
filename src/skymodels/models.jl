export SkyModel

abstract type AbstractSkyModel end

struct SkyModel{F,P,G<:AbstractDomain, A<:FourierTransform, M} <: AbstractSkyModel
    f::F
    prior::P
    grid::G
    algorithm::A
    metadata::M
end

function SkyModel(f, prior, grid::AbstractRectiGrid; algorithm = NFFTAlg(), metadata=nothing)
    return SkyModel(f, prior, grid, algorithm, metadata)
end

function VLBISkyModels.FourierDualDomain(grid::AbstractRectiGrid, array::AbstractArrayConfiguration, alg::FourierTransform; executor=Serial())
    return FourierDualDomain(grid, domain(array; executor), alg)
end

struct ObservedSkyModel{F, G<:VLBISkyModels.AbstractFourierDualDomain, M} <: AbstractSkyModel
    f::F
    grid::G
    metadata::M
end

function domain(m::AbstractSkyModel)
    return getfield(m, :grid)
end
function ObservedSkyModel(m::SkyModel, arr::AbstractArrayConfiguration)
    return ObservedSkyModel(m.f, FourierDualDomain(m.grid, arr, m.algorithm), m.metadata)
end

function set_array(m::AbstractSkyModel, array::AbstractArrayConfiguration)
    return ObservedSkyModel(m, array), m.prior
end


function idealvisibilities(m::AbstractSkyModel, x)
    skym = skymodel(m, x.sky)
    return visibilitymap(skym, domain(m))
end

function skymodel(m::AbstractSkyModel, x)
    return m.f(x, m.metadata)
end

preallocate_image(x, ::AbstractSingleDomain, ::FourierTransform) = x
function preallocate_image(x::AbstractRectiGrid, uv::AbstractSingleDomain, alg::FourierTransform)
    return FourierDualDomain(x, uv, alg)
end
