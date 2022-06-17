export stretched, shifted, rotated, renormed


"""
    $(TYPEDEF)
Abstract type for image modifiers. These are some model wrappers
that can transform any model using simple Fourier transform properties.
To see the implemented modifier
"""
abstract type AbstractModifier{M<:AbstractModel} <: AbstractModel end

"""
    $(SIGNATURES)
Returns the base model from a modified model.
"""
basemodel(model::AbstractModifier) = model.model
basemodel(model::AbstractModel) = model

unmodified(model::AbstractModel) = basemodel(model)
unmodified(model::AbstractModifier) = unmodified(basemodel(model))


flux(m::AbstractModifier) = flux(m.model)

Base.@aggressive_constprop @inline visanalytic(::Type{<:AbstractModifier{M}}) where {M} = visanalytic(M)
Base.@aggressive_constprop @inline imanalytic(::Type{<:AbstractModifier{M}}) where {M} = imanalytic(M)

radialextent(m::AbstractModifier) = radialextent(basemodel(m))


# @inline function apply_uv_transform(m::AbstractModifier, u::Number, v::Number, scale::Number)
#     ut, vt = transform_uv(m, u, v)
#     scale *= scale_uv(m, u, v)
#     return apply_uv_transform(basemodel(m), ut, vt, scale)
# end

# @inline function apply_uv_transform(::AbstractModel, u::Number, v::Number, scale::Number)
#     return u, v, scale
# end

# function apply_uv_transform(m::AbstractModifier, u::AbstractVector, v::AbstractVector)
#     ut = similar(u)
#     vt = similar(v)
#     T = typeof(scale_uv(m, 0.0, 0.0))
#     scale = ones(T, size(ut))
#     @inbounds for i in eachindex(u,v)
#         up, vp, sp = apply_uv_transform(m, u[i], v[i], scale[i])
#         ut[i] = up
#         vt[i] = vp
#         scale[i] = sp
#     end
#     return ut, vt, scale
# end

#@inline function _visibilities(m::AbstractModifier, u::AbstractArray, v::AbstractArray, args...)
#    ut, vt, scales = apply_uv_transform(m, u, v)
#    scales.*visibilities(unmodified(m), ut, vt, args...)
#end



function modelimage(::NotAnalytic,
    model::AbstractModifier,
    image::ComradeBase.AbstractIntensityMap, alg::FFTAlg, executor=SequentialEx())

    @set model.model = modelimage(model.model, image, alg, executor)
end

function modelimage(::NotAnalytic,
    model::AbstractModifier,
    image::ComradeBase.AbstractIntensityMap, alg::NUFT,
    executor=SequentialEx())
    _modelimage(model, image, alg, executor)
end



@inline function visibility_point(m::AbstractModifier, u, v, args...)
    ut, vt = transform_uv(m, u, v)
    scale = scale_uv(m, u, v)
    scale*visibility(basemodel(m), ut, vt, args...)
end

@inline function ComradeBase.intensity_point(m::AbstractModifier, x, y)
    xt, yt = transform_image(m, x, y)
    scale = scale_image(m, x, y)
    scale*ComradeBase.intensity_point(basemodel(m), xt, yt)
end

"""
    $(TYPEDEF)
Shifts the model by `Δx` units in the x-direction and `Δy` units
in the y-direction.
"""
struct ShiftedModel{M<:AbstractModel,T} <: AbstractModifier{M}
    model::M
    Δx::T
    Δy::T
end

"""
    $(SIGNATURES)
Shifts the model `m` in the image domain by an amount `Δx,Δy`.
"""
shifted(model, Δx, Δy) = ShiftedModel(model, Δx, Δy)
# This is a simple overload to simplify the type system
shifted(model::ShiftedModel, Δx, Δy) = ShiftedModel(basemodel(model), Δx+model.Δx, Δy+model.Δy)
radialextent(model::ShiftedModel, Δx, Δy) = radialextent(model.model) + hypot(abs(Δx), abs(Δy))

@inline transform_image(model::ShiftedModel, x, y) = (x-model.Δx, y-model.Δy)
@inline transform_uv(model::ShiftedModel, u, v) = (u, v)

@inline scale_image(model::ShiftedModel, x, y) = 1.0
@inline scale_uv(model::ShiftedModel, u, v) = exp(2im*π*(u*model.Δx + v*model.Δy))

@inline visanalytic(::Type{<:ShiftedModel{M}}) where {M} = visanalytic(M)


"""
    $(TYPEDEF)
Renormalizes the flux of the model to the new value `flux`.
We have also overloaded the Base.:* operator as syntactic sugar
although I may get rid of this.
"""
struct RenormalizedModel{M<:AbstractModel,T} <: AbstractModifier{M}
    model::M
    scale::T
    RenormalizedModel(model::M, f::T) where {M,T} = new{M,T}(model, f)
end

"""
    $(SIGNATURES)
Renormalizes the model `m` to have total flux `flux`.
"""
renormed(model::M, f) where {M<:AbstractModel} = RenormalizedModel(model, f)
Base.:*(model::AbstractModel, f::Real) = renormed(model, f)
Base.:*(f::Real, model::AbstractModel) = renormed(model, f)
Base.:/(f::Real, model::AbstractModel) = renormed(model, inv(f))
Base.:/(model::AbstractModel, f::Real) = renormed(model, inv(f))
# Dispatch on RenormalizedModel so that I just make a new RenormalizedModel with a different f
# This will make it easier on the compiler.
Base.:*(model::RenormalizedModel, f::Real) = renormed(model.model, model.scale*f)
# Overload the unary negation operator to be the same model with negative flux
Base.:-(model::AbstractModel) = renormed(model, -1.0)
flux(m::RenormalizedModel) = m.scale*flux(m.model)

@inline transform_image(::RenormalizedModel, x, y) = (x, y)
@inline transform_uv(::RenormalizedModel, u, v) = (u, v)

@inline scale_image(model::RenormalizedModel, x, y) = model.scale
@inline scale_uv(model::RenormalizedModel, u, v) = model.scale

#function _visibilities(m::RenormalizedModel, u::AbstractArray, v::AbstractArray)
#    m.scale*_visibilities(basemodel(m), u, v)
#end

@inline visanalytic(::Type{<:RenormalizedModel{M}}) where {M} = visanalytic(M)


"""
    $(TYPEDEF)
Stretched the model in the x and y directions, i.e. the new intensity is
```math
    I_s(x,y) = 1/(αβ) I(x/α, y/β),
```
where were renormalize the intensity to preserve the models flux.
"""
struct StretchedModel{M<:AbstractModel,T} <: AbstractModifier{M}
    model::M
    α::T
    β::T
end

"""
    $(SIGNATURES)
Stretches the model `m` according to the formula
```math
    I_s(x,y) = 1/(αβ) I(x/α, y/β),
```
where were renormalize the intensity to preserve the models flux.
"""
stretched(model, α, β) = StretchedModel(model, α, β)
radialextent(model::StretchedModel) = hypot(model.α, model.β)*radialextent(basemodel(model))

@inline transform_image(model::StretchedModel, x, y) = (x/model.α, y/model.β)
@inline transform_uv(model::StretchedModel, u, v) = (u*model.α, v*model.β)

@inline scale_image(model::StretchedModel, x, y) = inv(model.α*model.β)
@inline scale_uv(::StretchedModel, u, v) = one(eltype(u))

@inline visanalytic(::Type{<:StretchedModel{M}}) where {M} = visanalytic(M)


"""
    $(TYPEDEF)
Type for the rotated model. This is more fine grained constrol of
rotated model. In most use cases the end-user should be using
the `rotate` method e.g.

```julia
rotate(model, ξ)
```
"""
struct RotatedModel{M<:AbstractModel,T} <: AbstractModifier{M}
    model::M
    s::T
    c::T
end
function RotatedModel(model::T, ξ::F) where {T, F}
    s,c = sincos(ξ)
    return RotatedModel(model, s, c)
end

@inline visanalytic(::Type{<:RotatedModel{M}}) where {M} = visanalytic(M)

"""
    $(SIGNATURES)
Rotates the model by an amount `ξ` in radians.
"""
rotated(model, ξ) = RotatedModel(model, ξ)
posangle(model::RotatedModel) = atan(model.s, model.c)

@inline function transform_image(model::RotatedModel, x, y)
    s,c = model.s, model.c
    return c*x + s*y, -s*x + c*y
end

@inline function transform_uv(model::RotatedModel, u, v)
    s,c = model.s, model.c
    return c*u + s*v, -s*u + c*v
end

@inline scale_image(model::RotatedModel, x, y) = 1.0
@inline scale_uv(model::RotatedModel, u, v) = 1.0
