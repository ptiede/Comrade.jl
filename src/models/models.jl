export visibility, amplitude, closure_phase, logclosure_amplitude, bispectrum,
       visibilities, amplitudes, closure_phases, logclosure_amplitudes, bispectra,
       flux, intensitymap, intensitymap!, PolarizedModel
abstract type AbstractModel end
abstract type AbstactPolarizedModel <: AbstractModel end



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

include(joinpath(@__DIR__, "fft_alg.jl"))
include(joinpath(@__DIR__, "modelimage.jl"))
include(joinpath(@__DIR__, "modifiers.jl"))
include(joinpath(@__DIR__, "combinators.jl"))
include(joinpath(@__DIR__, "geometric_models.jl"))
include(joinpath(@__DIR__, "radio_image_models.jl"))








"""
    $(SIGNATURES)
Computes the complex visibility of model `m` at u,v positions `u,v`

If you want to compute the visibilities at a large number of positions
consider using the `visibilities` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function visibility(mimg::M, u, v) where {M}
    #first we split based on whether the model is primitive
    _visibility(isprimitive(M), mimg, u, v)
end

"""
    $(SIGNATURES)
Computes the visibility amplitude of model `m` at u,v positions `u,v`

If you want to compute the amplitudes at a large number of positions
consider using the `amplitudes` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function amplitude(model, u, v)
    return abs(visibility(model, u, v))
end

"""
    $(SIGNATURES)
Computes the complex bispectrum of model `m` at the uv-triangle
u1,v1 -> u2,v2 -> u3,v3

If you want to compute the bispectrum over a number of triangles
consider using the `bispectra` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function bispectrum(model, u1, v1, u2, v2, u3, v3)
    return visibility(model, u1, v1)*visibility(model, u2, v2)*visibility(model, u3, v3)
end

"""
    $(SIGNATURES)
Computes the closure phase of model `m` at the uv-triangle
u1,v1 -> u2,v2 -> u3,v3

If you want to compute closure phases over a number of triangles
consider using the `closure_phases` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function closure_phase(model, u1, v1, u2, v2, u3, v3)
    return angle(bispectrum(model, u1, v1, u2, v2, u3, v3))
end

"""
    $(SIGNATURES)
Computes the log-closure amplitude of model `m` at the uv-quadrangle
u1,v1 -> u2,v2 -> u3,v3 -> u4,v3 using the formula

```math
C = \\log\\left|\\frac{V(u1,v1)V(u2,v2)}{V(u3,v3)V(u4,v4)}\\right|
```

If you want to compute log closure amplitudes over a number of triangles
consider using the `logclosure_amplitudes` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function logclosure_amplitude(model, u1, v1, u2, v2, u3, v3, u4, v4)
    a1 = amplitude(model, u1, v1)
    a2 = amplitude(model, u2, v2)
    a3 = amplitude(model, u3, v3)
    a4 = amplitude(model, u4, v4)

    return log(a1*a2/(a3*a4))
end


#=
    Welcome to the trait jungle. What is below is
    how we specify how to evaluate the model
=#
@inline function _visibility(::NotPrimitive, m, u, v)
    return visibility_point(m, u, v)
end

@inline function _visibility(::IsPrimitive, m::M, u, v) where {M}
    _visibility_primitive(visanalytic(M), m, u, v)
end


@inline function _visibility_primitive(::IsAnalytic, mimg, u, v)
    return visibility_point(mimg, u, v)
end

@inline function _visibility_primitive(::NotAnalytic, mimg, u, v)
    return mimg.cache.sitp(u, v)
end



function visibilities(m, u::AbstractArray, v::AbstractArray)
    f(x,y) = visibility(m, x, y)
    return mappedarray(f, u, v)
end

function amplitudes(m, u::AbstractArray, v::AbstractArray)
    f(x,y) = amplitude(m, x, y)
    return mappedarray(f, u, v)
end

function bispectra(m,
                    u1::AbstractArray, v1::AbstractArray,
                    u2::AbstractArray, v2::AbstractArray,
                    u3::AbstractArray, v3::AbstractArray
                   )
    f(x1,y1,x2,y2,x3,y3) = bispectra(m, x1, y1, x2, y2, x3, y3)
    return mappedarray(f, u1, v1, u2, v2, u3, v3)
end

function closure_phases(m,
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray
                       )
    f(x1,y1,x2,y2,x3,y3) = closure_phase(m, x1, y1, x2, y2, x3, y3)
    return mappedarray(f, u1, v1, u2, v2, u3, v3)
end

function logclosure_amplitudes(m,
                               u1::AbstractArray, v1::AbstractArray,
                               u2::AbstractArray, v2::AbstractArray,
                               u3::AbstractArray, v3::AbstractArray,
                               u4::AbstractArray, v4::AbstractArray
                              )
    f(x1,y1,x2,y2,x3,y3,x4,y4) = logclosure_amplitude(m, x1, y1, x2, y2, x3, y3, x4, y4)
    return mappedarray(f, u1, v1, u2, v2, u3, v3, u4, v4)
end

function intensitymap!(im::IntensityMap, m::M) where {M}
    return intensitymap!(imanalytic(M), im, m)
end


function intensitymap!(::IsAnalytic, im::IntensityMap, m)
    xitr, yitr = imagepixels(im)
    @inbounds for (i,x) in pairs(xitr), (j,y) in pairs(yitr)
        im[j, i] = intensity_point(m, x, y)
    end
    return im
end


function intensitymap!(::NotAnalytic, img::IntensityMap, m)
    ny, nx = size(img)
    vis = ifftshift(phasedecenter!(fouriermap(m, img), img))
    ifft!(vis)
    for I in CartesianIndices(img)
        img[I] = real(vis[I])/(nx*ny)
    end
end
