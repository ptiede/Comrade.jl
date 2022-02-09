export modelimage

"""
    $(TYPEDEF)
Container for non-analytic model that contains a image cache which will hold the image,
a Fourier transform cache C, which usually an instance of a <: FourierCache which usually
holds a interpolator that allows you to compute visibilities.

# Notes
This is an internal implementation detail that shouldn't usually be called directly.
Instead the user should use the exported function `modelimage`, for example

```julia
using Comrade
m = ExtendedRing(20.0, 5.0)

# This creates an version where the image is dynamically specified according to the
# radial extent of the image
mimg = modelimage(m) # you can also optionally pass the number of pixels nx and ny

# Or you can create an IntensityMap
img = intensitymap(m, 100.0, 100.0, 512, 512)
mimg = modelimage(m, img)
```
"""
struct ModelImage{M,I,C} <: AbstractModelImage{M}
    model::M
    image::I
    cache::C
end
@inline visanalytic(::Type{<:ModelImage{M}}) where {M} = visanalytic(M)
@inline imanalytic(::Type{<:ModelImage{M}}) where {M} = imanalytic(M)
@inline isprimitive(::Type{<:ModelImage{M}}) where {M} = isprimitive(M)

function Base.show(io::IO, mi::ModelImage)
   #io = IOContext(io, :compact=>true)
   #s = summary(mi)
   ci = first(split(summary(mi.cache), "{"))
   println(io, "ModelImage")
   println(io, "\tmodel: ", summary(mi.model))
   println(io, "\timage: ", summary(mi.image))
   println(io, "\tcache: ", ci)
end

model(m::AbstractModelImage) = m.model
flux(mimg::ModelImage) = flux(mimg.image)

# function intensitymap(mimg::ModelImage)
#     intensitymap!(mimg.image, mimg.model)
#     mimg.image
# end

radialextent(m::ModelImage) = fov(m.image)[1]/2

#@inline visibility_point(m::AbstractModelImage, u, v) = visibility_point(model(m), u, v)

@inline intensity_point(m::AbstractModelImage, x, y) = intensity_point(model(m), x, y)

"""
    $(SIGNATURES)
Construct a `ModelImage` from a `model`, `image` and the optionally
specified visibility algorithm `alg`

# Notes
For analytic models this is a no-op and just return the model.
For non-analytic models this wraps the model in a object with an image
and precomputes the fourier transform using `alg`.
"""
@inline function modelimage(model::M, image::ComradeBase.AbstractIntensityMap, alg::FourierTransform=FFTAlg()) where {M}
    return modelimage(visanalytic(M), model, image, alg)
end

@inline function modelimage(::IsAnalytic, model, image, args...; kwargs...)
    return model
end

function _modelimage(model, image, alg)
    intensitymap!(image, model)
    cache = create_cache(alg, image)
    return ModelImage(model, image, cache)
end

@inline function modelimage(::NotAnalytic, model, image::ComradeBase.AbstractIntensityMap, alg::FourierTransform=FFTAlg())
    _modelimage(model, image, alg)
end

@inline function modelimage(model::M, cache::AbstractCache) where {M}
    return modelimage(visanalytic(M), model, cache)
end

@inline function modelimage(::IsAnalytic, model, cache::AbstractCache)
    return model
end

@inline function modelimage(::NotAnalytic, model, cache::AbstractCache)
    img = cache.pimg
    intensitymap!(img, model)
    newcache = updatecache(cache, img)
    return ModelImage(model, img, newcache)
end

"""
    $(SIGNATURES)
Construct a `ModelImage` where just the model `m` is specified

# Notes
If m `IsAnalytic()` is the visibility domain this is a no-op and just returns the model itself.
Otherwise `modelimage` will *guess* a reasonable field of view based on the `radialextent`
function. One can optionally pass the number of pixels nx and ny in each direction.
"""
function modelimage(m::M;
                    fovx=2*radialextent(m),
                    fovy=2*radialextent(m),
                    nx=512,
                    ny=512,
                    pulse=ComradeBase.DeltaPulse(),
                    alg=FFTAlg()) where {M}
    if visanalytic(M) == IsAnalytic()
        return m
    else
        T = typeof(intensity_point(m, 0.0, 0.0))
        img = IntensityMap(zeros(T,ny,nx), fovx, fovy, pulse)
        modelimage(m, img, alg)
    end
end
