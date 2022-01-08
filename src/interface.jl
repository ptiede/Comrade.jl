abstract type AbstractModel end

"""
    $(TYPEDEF)
This trait specifies whether the model is a *primitive*

# Notes
This will likely turn into a trait in the future so people
can inject their models into ROSEx more easily.
"""
abstract type PrimitiveTrait end
struct IsPrimitive end
struct NotPrimitive end

"""
    isprimitive(::Type)
Dispatch function that specifies whether a type is a primitive ROSEx model.
This function is used for dispatch purposes when composing models.

# Notes
If a user is specifying their own model primitive model outside of ROSEx they need
to specify if it is primitive

```julia
struct MyPrimitiveModel end
ROSEx.isprimitive(::Type{MyModel}) = ROSEx.IsPrimitive()
```

"""
function isprimitive end

isprimitive(::Type{<:AbstractModel}) = NotPrimitive()

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
Computes the complex visibility of model `m` at u,v positions `u,v`

If you want to compute the visibilities at a large number of positions
consider using the `visibilities` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
function visibility end
