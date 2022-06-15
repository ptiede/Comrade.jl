abstract type AbstractModel end
abstract type AbstractPolarizedModel <: AbstractModel end




"""
    $(TYPEDEF)
This trait specifies whether the model is a *primitive*

# Notes
This will likely turn into a trait in the future so people
can inject their models into Comrade more easily.
"""
abstract type PrimitiveTrait end

"""
    $(TYPEDEF)
Trait for primitive model
"""
struct IsPrimitive end
"""
    $(TYPEDEF)
Trait for not-primitive model
"""
struct NotPrimitive end

"""
    isprimitive(::Type)
Dispatch function that specifies whether a type is a primitive Comrade model.
This function is used for dispatch purposes when composing models.

# Notes
If a user is specifying their own model primitive model outside of Comrade they need
to specify if it is primitive

```julia
struct MyPrimitiveModel end
ComradeBase.isprimitive(::Type{MyModel}) = ComradeBase.IsPrimitive()
```

"""
function isprimitive end

@inline isprimitive(::Type{<:AbstractModel}) = NotPrimitive()

"""
    DensityAnalytic
Internal type for specifying the nature of the model functions.
Whether they can be easily evaluated pointwise analytic. This
is an internal type that may change.
"""
abstract type DensityAnalytic end

"""
    $(TYPEDEF)
Defines a trait that a states that a model is analytic.
This is usually used with an abstract model where we use
it to specify whether a model has a analytic fourier transform
and/or image.
"""
struct IsAnalytic <: DensityAnalytic end

"""
    $(TYPEDEF)
Defines a trait that a states that a model is analytic.
This is usually used with an abstract model where we use
it to specify whether a model has does not have a easy analytic
fourier transform and/or intensity function.
"""
struct NotAnalytic <: DensityAnalytic end

"""
    visanalytic(::Type{<:AbstractModel})
Determines whether the model is pointwise analytic in Fourier domain, i.e. we can evaluate
its fourier transform at an arbritrary point.

If `IsAnalytic()` then it will try to call `visibility_point` to calculate the complex visibilities.
Otherwise it fallback to using the FFT that works for all models that can compute an image.

"""
@inline visanalytic(::Type{<:AbstractModel}) = NotAnalytic()

"""
    imanalytic(::Type{<:AbstractModel})
Determines whether the model is pointwise analytic in the image domain, i.e. we can evaluate
its intensity at an arbritrary point.

If `IsAnalytic()` then it will try to call `intensity_point` to calculate the intensity.
"""
@inline imanalytic(::Type{<:AbstractModel}) = IsAnalytic()



#=
    This is internal function definitions for how to
    compose whether a model is analytic. We need this
    for composite models.
=#
@inline Base.:*(::IsAnalytic, ::IsAnalytic) = IsAnalytic()
@inline Base.:*(::IsAnalytic, ::NotAnalytic) = NotAnalytic()
@inline Base.:*(::NotAnalytic, ::IsAnalytic) = NotAnalytic()
@inline Base.:*(::NotAnalytic, ::NotAnalytic) = NotAnalytic()


"""
    $(SIGNATURES)
Function that computes the pointwise visibility if the model has the trait
in the fourier domain `IsAnalytic()`. Otherwise it will use the FFTW fallback.
"""
function visibility_point end

"""
    $(SIGNATURES)
Function that computes the pointwise intensity if the model has the trait
in the image domain `IsAnalytic()`. Otherwise it will use construct the image in visibility
space and invert it.
"""
function intensity_point end


"""
    $(SIGNATURES)
Computes the intensity map of model by modifying the input IntensityMap object
"""
function intensitymap! end


"""
    $(SIGNATURES)
Computes the intensity map of model. This version requires additional information to
construct the grid.

# Example

```julia
m = Gaussian()
# field of view
fovx, fovy = 5.0
fovy = 5.0
# number of pixels
nx, ny = 128

img = intensitymap(m, fovx, fovy, nx, ny; pulse=DeltaPulse())
```
"""
function intensitymap end


# """
#     ($SIGNATURES)
# A threaded version of the intensitymap function.

# # Notes
# If using autodiff this won't play well with zygote given the mutation. This will be fixed in
# future versions.
# """
# function tintensitymap end

# """
#     ($SIGNATURES)
# A threaded version of the intensitymap! function.

# # Notes
# If using autodiff this won't play well with zygote given the mutation. This will be fixed in
# future versions.
# """
# function tintensitymap! end
