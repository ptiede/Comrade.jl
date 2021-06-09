"""
    $(TYPEDEF)
Abstract type for image modifiers. These are some model wrappers
that can transform any model using simple Fourier transform properties.
To see the implemented modifier
"""
abstract type AbstractModifier{M<:AbstractModel,T} <: AbstractModel{T} end

ImStyle(::Type{<:AbstractModifier{M,T}}) where {M,T} = ImStyle(M)
VisStyle(::Type{<:AbstractModifier{M,T}}) where {M,T} = VisAnalytic()

"""
    $(SIGNATURES)
Returns a list containing all the image modifiers available in ROSE.
"""
modifierlist() = subtypes(AbstractModifier)

VisStyle(::Type{<:AbstractModifier}) = VisAnalytic()

@inline transform(m::AbstractModifier, x, y) = (x,y)
@inline imscale(m::AbstractModifier, x, y) = one(eltype(x))

@inline transformuv(m::AbstractModifier, u, v)  = (u,v)
@inline visscale(m::AbstractModifier, u, v) = one(eltype(u))

@inline function intensity(::ImAnalytic, m::AbstractModifier, x, y, args...)
    return imscale(m,x,y)*intensity(basemodel(m), transform(m, x, y), args...)
end

@inline function visibility(m::AbstractModifier, u, v, args...)
    return visscale(m, u, v)*visibility(basemodel(m), transformuv(m, u, v), args...)
end


# TODO: Make basemodel a trait!
# TODO: How should I deal with nested modifiers?
"""
    $(SIGNATURES)
Returns the base model from a modified model.
"""
basemodel(model::AbstractModifier) = model.model
flux(m::AbstractModifier) = flux(m.model)




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
struct SmoothedModel{M<:AbstractModel, T} <: AbstractModel{T}
    model::M
    σ::T
end
basemodel(m::SmoothedModel) = m.model
smoothed(model, σ) = SmoothedModel(model, σ)
ImStyle(::Type{SmoothedModel{M,T}}) where {M,T} = ImNumeric()

visscale(m::SmoothedModel, u, v) = exp(-2*(π*m.σ)^2*(u^2+v^2))

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

function visibilities(::VisAnalytic, m::SmoothedModel, cache, args...)
    vis = visibilities(basemodel(m), cache, args...)
    u, v = cache.u, cache.v
    for i in eachindex(u,v)
        vis[i] *= visscale(m, u[i], v[i])
    end
    return vis
end

function visibilities(::VisAnalytic, m::SmoothedModel, u, v, args...)
    vis = visibilities(basemodel(m), u, v, args...)
    for i in eachindex(u,v)
        vis[i] *= visscale(m, u[i], v[i])
    end
    return vis
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

@inline transform(m::ShiftedModel, x, y) = x-m.Δx, y-m.Δy
@inline visscale(m::ShiftedModel, u, v) = exp(2im*π*(u*m.Δx + v*m.Δy))



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

@inline visscale(m::RenormalizedModel, u, v) = m.flux/flux(basemodel(m))
@inline imscale(m::RenormalizedModel, x,y) = m.flux/flux(basemodel(m))



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

@inline transformuv(m::StretchedModel, u, v) = u*m.α, v*m.β
@inline transform(m::StretchedModel, x, y) = x/m.α, y/m.β
@inline imscale(m::StretchedModel, x, y) = 1/(m.α*m.β)



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

function transformuv(m::RotatedModel, u, v)
    s,c = m.s, m.c
    return c*u - s*v, s*u + c*v
end

function transform(m::RotatedModel, x, y)
    s,c = m.s, m.c
    return c*x - s*y, s*x + c*y
end
