


# TODO: Make basemodel a trait!
# TODO: How should I deal with nested modifiers?
"""
    $(SIGNATURES)
Returns the base model from a modified model.
"""
basemodel(model::AbstractModifier) = model.model
flux(m::AbstractModifier) = flux(m.model)


abstract type CompositeModel{T} <: AbstractModel{T} end

"""
    $(TYPEDEF)
Adds two models together to create composite models. Note that
I may change this in the future so make it easier on the compiler,
i.e. make the composite model a fancy model vector and heap allocate
stuff. This should help when combining multiple models together.
"""
struct AddModel{T,T1<:AbstractModel{T},T2<:AbstractModel{T}} <: CompositeModel{T}
    m1::T1
    m2::T2
end
Base.:+(m1::T1, m2::T2) where {T1<:AbstractModel, T2<:AbstractModel} = AddModel(m1, m2)
add(m1::M1, m2::M2) where {M1<:AbstractModel, M2<:AbstractModel} = AddModel(m1, m2)

"""
    $(SIGNATURES)
Returns the components for a composite model. This
will return a Tuple with all the models you have constructed.
"""
components(m::AbstractModel) = m
components(m::AddModel{M1,M2}) where
    {M1<:AbstractModel, M2<:AbstractModel} = (components(m.m1)..., components(m.m2)...)

flux(m::AddModel) = flux(m.m1) + flux(m.m2)

function intensity(m::AddModel{T1,T2}, x, y, args...) where {T1, T2}
    return intensity(m.m1, x, y, args...) + intensity(m.m2, x, y, args...)
end


@inline function visibility(m::AddModel{T1,T2}, u, v, args...) where {T1, T2}

    return visibility(m.m1, u, v, args...) + visibility(m.m2, u, v, args...)
end

"""
    $(TYPEDEF)
Smooths a `model` with a Gaussian kernel with standard deviation
of `σ`.

## Notes

```julia
m = smoothed(Disk(), 5.0)
intensity(m, 2.0, 2.0; fov=10.0, npix=128)
```
will compute the intensity of the smoothed disk with using an interpolation
with `128` pixel nodes and a total `fov` of 10 in x and y direction.

**This needs to be improved**
"""
struct SmoothedModel{M<:AbstractModel, T} <: AbstractModifier{M,T}
    model::M
    σ::T
end
smoothed(model, σ) = SmoothedModel(model, σ)
ImStyle(::Type{SmoothedModel{M,T}}) where {M,T} = ImGrid()


function intensitymap!(sim::StokesImage{T,S}, model::SmoothedModel) where {T,S}
    cim = deepcopy(sim)
    intensitymap!(cim, basemodel(model))
    dx,dy = pixelsizes(sim)
    xitr,yitr = imagepixels(sim)
    σ_px = model.σ/abs(dx)

    # Now I need to pick my kernel size. I am going out to 5σ for the
    # gaussian kernel. I have to add one for the convolution to play nice
    nkern = Int(floor(σ_px)*10 + 1)

    imfilter!(sim, cim,
        gaussian((σ_px, σ_px),(nkern,nkern)),
        Fill(0.0, cim),
        FFT())
    return sim
end


@inline function visibility(model::SmoothedModel, u, v, args...)
    return visibility(basemodel(model), u, v, args...)*
            exp(-2*(π*model.σ)^2*(u^2+v^2))
end

"""
    $(TYPEDEF)
Shifts the model by `Δx` units in the x-direction and `Δy` units
in the y-direction.
"""
struct ShiftedModel{T,M<:AbstractModel} <: AbstractModifier{M,T}
    model::M
    Δx::T
    Δy::T
end
shifted(model, Δx, Δy) = ShiftedModel(model, Δx, Δy)

@inline function intensity(model::ShiftedModel, x, y, args...)
    return intensity(basemodel(model), x-model.Δx, y-model.Δy, args...)
end

@inline function visibility(model::ShiftedModel, u, v, args...)
    return visibility(basemodel(model), u, v, args...)*
            exp(2im*π*(u*model.Δx + v*model.Δy))
end

"""
    $(TYPEDEF)
Renormalizes the flux of the model to the new value `flux`.
We have also overloaded the Base.:* operator as syntactic sugar
although I may get rid of this.
"""
struct RenormalizedModel{M<:AbstractModel,T} <: AbstractModifier{M,T}
    model::M
    flux::T
end
renormed(model::M, flux) where {M<:AbstractModel} = RenormalizedModel(model, flux)
Base.:*(model::AbstractModel, flux::Real) = renormed(model, flux)
# Dispatch on RenormalizedModel so that I just make a new RenormalizedModel with a different flux
# This will make it easier on the compiler.
Base.:*(model::RenormalizedModel, flux::Real) = renormed(model.model, model.flux*flux)
Base.:*(flux::Real, model::AbstractModel) = renormed(model, flux)
# Overload the unary negation operator to be the same model with negative flux
Base.:-(model::AbstractModel) = renormed(model, -flux(model))
flux(m::RenormalizedModel) = m.flux


@inline function intensity(m::RenormalizedModel, x, y, args...)
    return intensity(basemodel(m), x,y, args...)*m.flux/flux(basemodel(m))
end

@inline function visibility(model::RenormalizedModel,
                            u, v, args...)
    return visibility(basemodel(model), u, v, args...)*model.flux/flux(basemodel(model))
end

"""
    $(TYPEDEF)
Stretched the model in the x and y directions, i.e. the new intensity is
```math
    I_s(x,y) = 1/(αβ) I(x/α, y/β),
```
where were renormalize the intensity to preserve the models flux.
"""
struct StretchedModel{M<:AbstractModel,T} <: AbstractModifier{M,T}
    model::M
    α::T
    β::T
end
stretched(model, α, β) = StretchedModel(model, α, β)

@inline function intensity(model::StretchedModel, x,y, args...)
    return intensity(basemodel(model), x/model.α, y/model.β, args...)/(model.α*model.β)
end
@inline function visibility(model::StretchedModel, u, v, args...)
    return visibility(basemodel(model), u*model.α, v*model.β, args...)
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
struct RotatedModel{M<:AbstractModel,T} <: AbstractModifier{M,T}
    model::M
    s::T
    c::T
end
function RotatedModel(model::T, ξ::F) where {T, F}
    s,c = sincos(ξ)
    return RotatedModel(model, s, c)
end
rotated(model, ξ) = RotatedModel(model, ξ)
posangle(model::RotatedModel) = atan(model.s, model.c)


@inline function intensity(model::RotatedModel, x,y, args...)
    s,c = model.s, model.c
    xx, yy = c*x - s*y, s*x + c*y
    return intensity(basemodel(model), xx, yy, args...)
end

@inline function visibility(model::RotatedModel, u,v, args...)
    s,c = model.s, model.c
    uu, vv = c*u - s*v, s*u + c*v
    return visibility(basemodel(model), uu, vv, args...)
end
