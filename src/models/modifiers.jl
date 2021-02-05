


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
Adds two models together to create composite models. Note that
I may change this in the future so make it easier on the compiler,
i.e. make the composite model a fancy model vector and heap allocate
stuff. This should help when combining multiple models together.
"""
struct AddModel{T,T1<:AbstractModel{T},T2<:AbstractModel{T}} <: AbstractModifier{T}
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

function visibility(::IsAnalytic, m::AddModel{T1,T2}, u, v, args...) where {T1, T2}
    return visibility(m.m1, u, v, args...) + visibility(m.m2, u, v, args...)
end

"""
    $(TYPEDEF)
Smooths a `model` with a Gaussian kernel with standard deviation
of `σ`.

## Notes
While the visibility representation of smoothing is simple due to the
convolution theorem, computing the intensities is not. To create
the smoothed model intensities I thus, use ImageFiltering to smooth
the pixels and then create a interpolator holding the results. This
interpolator is then cached using Memoization.jl to prevent the
computation from being repeated. As a result I need to first create
a image raster with some FOV and resolution. To do this when calling
intensity you can specify the fov and number of pixels

```julia
m = smoothed(Disk(), 5.0)
intensity(m, 2.0, 2.0; fov=10.0, npix=128)
```
will compute the intensity of the smoothed disk with using an interpolation
with `128` pixel nodes and a total `fov` of 10 in x and y direction.

**This needs to be improved**
"""
struct SmoothedModel{T<:AbstractModel, F} <: AbstractModifier{F}
    model::T
    σ::F
end
smoothed(model, σ) = SmoothedModel(model, σ)

Memoization.@memoize function _smoothitp(sim::StokesMatrix, σ)
    ny,nx = size(sim)
    dx,dy = pixelsizes(sim)
    xitr,yitr = pixel_iterator(sim)
    σ_px = σ/abs(dx)

    # Now I need to pick my kernel size. I am going out to 5σ for the
    # gaussian kernel. I have to add one for the convolution to play nice
    nkern = Int(floor(σ_px)*10 + 1)

    bsim = imfilter(sim,
        gaussian((σ_px, σ_px),(nkern,nkern)),
        Fill(0.0, sim),
        FFT())

    itp = interpolate(bsim, BSpline(Cubic(Line(OnGrid()))))
	etp = extrapolate(itp, 0)
    sitp = scale(etp, xitr, yitr)

    return sitp
end

Memoization.@memoize function _create_image(model, fov, npix)
    sim = stokesmatrix(model, npix, npix, fov, fov)
    return sim
end

function intensity(model::SmoothedModel, x, y, args...;fov=200.0, npix=512)
    sim = _create_image(basemodel(model), fov, npix)
    bitp = _smoothitp(sim, model.σ)
    return bitp(y, x)
end

function visibility(::IsAnalytic, model::SmoothedModel, u, v, args...)
    return visibility(basemodel(model), u, v, args...)*
            exp(-2*(π*model.σ)^2*(u^2+v^2))
end

"""
    $(TYPEDEF)
Shifts the model by `Δx` units in the x-direction and `Δy` units
in the y-direction.
"""
struct ShiftedModel{T<:AbstractModel,F} <: AbstractModifier{F}
    model::T
    Δx::F
    Δy::F
end
shifted(model, Δx, Δy) = ShiftedModel(model, Δx, Δy)

function intensity(model::ShiftedModel, x, y, args...)
    return intensity(basemodel(model), x-model.Δx, y-model.Δy, args...)
end

function visibility(::IsAnalytic, model::ShiftedModel, u, v, args...)
    return visibility(basemodel(model), u, v, args...)*
            exp(2im*π*(u*model.Δx + v*model.Δy))
end

"""
    $(TYPEDEF)
Renormalizes the flux of the model to the new value `flux`.
We have also overloaded the Base.:* operator as syntactic sugar
although I may get rid of this.
"""
struct RenormalizedModel{T<:AbstractModel,F} <: AbstractModifier{F}
    model::T
    flux::F
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


@inline function intensity(m::RenormalizedModel{T,F}, x, y, args...) where {T<:AbstractModel, F}
    return intensity(basemodel(m), x,y, args...)*m.flux/flux(basemodel(m))
end

@inline function visibility(::IsAnalytic,
                            model::RenormalizedModel{T,F},
                            u, v, args...) where {T<:AbstractModel, F}
    return visibility(basemodel(model), u, v, args...)*model.flux/flux(basemodel(model))
end

"""
    $(TYPEDEF)
Stretched the model in the x and y directions, i.e. the new intensity is
```math
    I_s(x,y) = αβI(αx, βy),
```
where were renormalize the intensity to preserve the models flux.
"""
struct StretchedModel{T<:AbstractModel,F} <: AbstractModifier{F}
    model::T
    α::F
    β::F
end
stretched(model, α, β) = StretchedModel(model, 1/α, 1/β)

function intensity(model::StretchedModel, x,y, args...)
    return intensity(basemodel(model), x*model.α, y*model.β, args...)*model.α*model.β
end
function visibility(::IsAnalytic, model::StretchedModel, u, v, args...)
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
struct RotatedModel{T<:AbstractModel,F} <: AbstractModifier{F}
    model::T
    s::F
    c::F
end
function RotatedModel(model::T, ξ::F) where {T, F}
    s,c = sincos(ξ)
    return RotatedModel(model, s, c)
end
rotated(model, ξ) = RotatedModel(model, ξ)
angle(model::RotatedModel) = atan(model.s, model.s)


function intensity(model::RotatedModel, x,y, args...)
    s,c = model.s, model.c
    xx, yy = c*x - s*y, s*x + c*y
    return intensity(basemodel(model), xx, yy, args...)
end

function visibility(::IsAnalytic ,model::RotatedModel, u,v, args...)
    s,c = model.s, model.c
    uu, vv = c*u - s*v, s*u + c*v
    return visibility(basemodel(model), uu, vv, args...)
end
