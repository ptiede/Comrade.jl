using LoopVectorization

abstract type AbstractModel end
abstract type AbstractRadioImage <: AbstractModel end
abstract type AbstractModifier <: AbstractModel end


struct Gaussian <: AbstractModel end
struct Disk{T} <: AbstractModel
    radius::T
    shift::T
end

struct Themage{S,M<:AbstractMatrix{S}} <: AbstractRadioImage
    coeff::M
    smooth::S
    psizex::S
    psizey::S
    function Themage(coeff::M, smooth::S) where {M,S}
        ny, nx = size(coeff)
        psizex = one(S)/max(nx-1, 1)
        psizey = one(S)/max(ny-1, 1)
        new{S,M}(coeff, smooth, psizex, psizey)
    end

end

@inline Base.size(model::Themage) = size(model.coeff)

function _basis(model::Themage, x::T) where {T}
    return ifelse(abs(x) < model.smooth, sinc(x)*sinc(x/model.smooth), zero(T))
    #return exp(-x^2/(2*model.smooth^2))
end

function intensity(model::Themage, x, y)
    sum = zero(eltype(model.coeff))
    ny,nx = size(model)
    @avx for I in CartesianIndices(model.coeff)
        dx = x*nx - (-nx/2 + I[2])
        dy = y*ny - (-ny/2 + I[1])
        sum += model.coeff[I]*_basis(model, dx)*_basis(model, dy)
    end
    return sum
end

function test(model::Themage)
    sum = 0.0
    @avx for i in eachindex(model.coeff)
        sum += model.coeff[i]
    end
    return sum
end

"""
    Returns the intensity for the model with args...
"""
function intensity(model::AbstractModel, args...) end

"""
    Returns the visibility for the model with args...
"""
function modelvisiblity(model::AbstractModel, args...) end


function intensity(model::Gaussian, x,y, args...)
    return exp(-sum(abs2,(x,y)))/sqrt(2π)
end

function visibility(model::Gaussian, u, v, args...)
    return exp(-2π^2*(u^2 + v^2))
end

basemodel(model::AbstractModifier) = model.model

struct ShiftModel{T<:AbstractModel,F} <: AbstractModifier
    model::T
    Δx::F
    Δy::F
end
shift(model, Δx, Δy) = ShiftModel(model, Δx, Δy)

struct ScaledModel{T<:AbstractModel,F} <: AbstractModifier
    model::T
    α::F
    β::F
end
scale(model, α, β) = ScaledModel(model, 1/α, 1/β)

struct RotateModel{T<:AbstractModel,F} <: AbstractModifier
    model::T
    s::F
    c::F
end
function RotateModel(model::T, ξ::F) where {T, F}
    s,c = sincos(ξ)
    return RotateModel(model, s, c)
end

function intensity(model::ShiftModel, x,y, args...)
    return intensity(basemodel(model), x-model.Δx, y-model.Δy, args...)
end
function visibility(model::ShiftModel, u, v, args...)
    return visibility(basemodel(model), u, v)*exp(-2π*im*(u*model.Δx + v*model.Δy))
end

function intensity(model::ScaledModel, x,y, args...)
    return intensity(basemodel(model), x*model.α, y*model.β, args...)*model.α*model.β
end
function visibility(model::ScaledModel, u, v, args...)
    return visibility(basemodel(model), u/model.α, v/model.β, args...)
end

function intensity(model::RotateModel, x,y, args...)
    s,c = model.s, model.c
    xx, yy = c*x - s*y, s*x + c*y
    return intensity(basemodel(model), xx, yy, args...)
end
function visibility(model::RotateModel, u,v, args...)
    s,c = model.s, model.c
    uu, vv = c*u - s*v, s*u + c*v
    return visibility(basemodel(model), uu, vv, args...)
end
