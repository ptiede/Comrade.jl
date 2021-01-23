basemodel(model::AbstractModifier) = model.model
flux(m::AbstractModifier) = flux(m.model)


struct AddModel{T,T1<:AbstractModel{T},T2<:AbstractModel{T}} <: AbstractModel{T}
    m1::T1
    m2::T2
end
Base.:+(m1::T1, m2::T2) where {T1<:AbstractModel, T2<:AbstractModel} = AddModel(m1, m2)

function intensity(m::AddModel{T1,T2}, x, y, args...) where {T1, T2}
    return intensity(m.m1, x, y, args...) + intensity(m.m1, x, y, args...)
end

function visibility(m::AddModel{T1,T2}, u, v, args...) where {T1, T2}
    return visibility(m.m1, u, v, args...) + visibilty(m.m1, u, v, args...)
end

struct ShiftModel{T<:AbstractModel,F} <: AbstractModifier{F}
    model::T
    Δx::F
    Δy::F
end
shift(model, Δx, Δy) = ShiftModel(model, Δx, Δy)

function intensity(model::ShiftModel, x, y, args...)
    return intensity(basemodel(model), x-model.Δx, y-model.Δy, args...)
end

function visibility(::IsAnalytic, model::ShiftModel, u, v)
    return visibility(basemodel(model), u, v, args...)*
            exp(-2im*π*(u*model.Δx + v*model.Δy))
end


struct ScaledModel{T<:AbstractModel,F} <: AbstractModifier{F}
    model::T
    α::F
    β::F
end
scale(model, α, β) = ScaledModel(model, 1/α, 1/β)

function intensity(model::ScaledModel, x,y, args...)
    return intensity(basemodel(model), x*model.α, y*model.β, args...)*model.α*model.β
end
function visibility(::IsAnalytic, model::ScaledModel, u, v, args...)
    return visibility(basemodel(model), u/model.α, v/model.β, args...)
end

"""
    $(TYPEDEF)
Type for the rotated model. This is more fine grained constrol of
rotated model. In most use cases the end-user should be using
the `rotate` method e.g.

```julia
rotate(model, ξ)
```
"""
struct RotateModel{T<:AbstractModel,F} <: AbstractModifier{F}
    model::T
    s::F
    c::F
end
function RotateModel(model::T, ξ::F) where {T, F}
    s,c = sincos(ξ)
    return RotateModel(model, s, c)
end
rotate(model, ξ) = RotateModel(model, ξ)
angle(model::RotateModel) = atan(m.s, m.s)


function intensity(model::RotateModel, x,y, args...)
    s,c = model.s, model.c
    xx, yy = c*x - s*y, s*x + c*y
    return intensity(basemodel(model), xx, yy, args...)
end

function visibility(::IsAnalytic ,model::RotateModel, u,v, args...)
    s,c = model.s, model.c
    uu, vv = c*u - s*v, s*u + c*v
    return visibility(basemodel(model), uu, vv, args...)
end
