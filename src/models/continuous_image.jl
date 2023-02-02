export ContinuousImage

"""
    ContinuousImage{A<:IntensityMap, P} <: AbstractModel
    ContinuousImage(img::Intensitymap, kernel)

The basic continuous image model for Comrade. This expects a IntensityMap style object as its imag
as well as a image kernel or pulse that allows you to evaluate the image at any image
and visibility location. The image model is

    I(x,y) = ∑ᵢ Iᵢⱼ κ(x-xᵢ, y-yᵢ)

where `Iᵢⱼ` are the flux densities of the image `img` and κ is the intensity function for the
`kernel`.


"""
struct ContinuousImage{A <: IntensityMapTypes, P} <: AbstractModel
    """
    Discrete representation of the image.
    """
    img::A
    """
    Image Kernel that transforms from the discrete image to a continuous one. This is
    sometimes called a pulse function in `eht-imaging`.
    """
    kernel::P
end

ComradeBase.ispolarized(::Type{<:ContinuousImage{A}}) where {A<:StokesIntensityMap} = IsPolarized()
ComradeBase.ispolarized(::Type{<:ContinuousImage{A}}) where {A<:IntensityMap{<:StokesParams}} = IsPolarized()
ComradeBase.ispolarized(::Type{<:ContinuousImage{A}}) where {A<:IntensityMap{<:Real}} = NotPolarized()

stokes(cimg::ContinuousImage, v) = ContinuousImage(stokes(parent(cimg), v), cimg.kernel)

Base.parent(m::ContinuousImage)         = m.img
Base.length(m::ContinuousImage)         = length(parent(m))
Base.size(m::ContinuousImage)           = size(parent(m))
Base.size(m::ContinuousImage, i::Int)   = size(parent(m), i::Int)
Base.firstindex(m::ContinuousImage)     = firstindex(parent(m))
Base.lastindex(m::ContinuousImage)      = lastindex(parent(m))
Base.iterate(m::ContinuousImage)        = iterate(parent(m))
Base.iterate(m::ContinuousImage, state) = iterate(parent(m), state)

Base.IteratorSize(::ContinuousImage{A, P}) where {A,P} = Base.IteratorSize(M)
Base.IteratorEltype(::ContinuousImage{A, P}) where {A,P} = Base.IteratorEltype(M)
Base.eltype(::ContinuousImage{A, P}) where {A,P} = eltype(A)

Base.getindex(img::ContinuousImage, args...) = getindex(parent(img), args...)
Base.axes(m::ContinuousImage) = axes(parent(m))
ComradeBase.imagegrid(m::ContinuousImage) = imagegrid(parent(m))
ComradeBase.AxisKeys.named_axiskeys(m::ContinuousImage) = named_axiskeys(parent(m))
ComradeBase.AxisKeys.axiskeys(m::ContinuousImage)       = ComradeBase.AxisKeys.axiskeys(parent(m))
ComradeBase.axisdims(m::ContinuousImage) = axisdims(parent(m))

Base.similar(m::ContinuousImage, ::Type{S}, dims) where {S} = ContinuousImage(similar(parent(m), S, dims), m.kernel)

function ContinuousImage(img::IntensityMapTypes, pulse::Pulse)
    return ContinuousImage{typeof(img), typeof(pulse)}(img, pulse)
end

ContinuousImage(img::IntensityMapTypes, cache::AbstractCache) = modelimage(ContinuousImage(img, cache.pulse), cache)

function ContinuousImage(im::AbstractMatrix, fovx::Real, fovy::Real, x0::Real, y0::Real, pulse, header=nothing)
    xitr, yitr = imagepixels(fovx, fovy, size(img, 1), size(img,2), x0, y0)
    img = IntensityMap(im, (X=xitr, Y=yitr), header)
    # spulse = stretched(pulse, step(xitr), step(yitr))
    return ContinuousImage(img, pulse)
end

function ContinuousImage(im::AbstractMatrix, fov::Real, x0::Real, y0::Real, pulse, header=nothing)
    return ContinuousImage(im, fov, fov, x0, y0, pulse, header)
end



ComradeBase.imagepixels(img::ContinuousImage) = NamedTuple{names.(dims(img.img))}(dims(img.img))

# IntensityMap will obey the Comrade interface. This is so I can make easy models
visanalytic(::Type{<:ContinuousImage}) = NotAnalytic() # not analytic b/c we want to hook into FFT stuff
imanalytic(::Type{<:ContinuousImage}) = IsAnalytic()
isprimitive(::Type{<:ContinuousImage}) = IsPrimitive()

radialextent(c::ContinuousImage) = maximum(values(fieldofview(c.img)))/2

function intensity_point(m::ContinuousImage, p)
    dx, dy = pixelsizes(m.img)
    sum = zero(eltype(m.img))
    ms = stretched(m.kernel, dx, dy)
    @inbounds for (I, p0) in pairs(imagegrid(m.img))
        dp = (X=(p.X - p0.X), Y=(p.Y - p0.Y))
        k = intensity_point(ms, dp)
        sum += m.img[I]*k
    end
    return sum
end

convolved(cimg::ContinuousImage, m::AbstractModel) = ContinuousImage(cimg.img, convolved(cimg.pulse, m))


"""
    modelimage(img::ContinuousImage, alg=NFFTAlg())

Create a model image directly using an image, i.e. treating it as the model. You
can optionally specify the Fourier transform algorithm using `alg`
"""
@inline function modelimage(model::ContinuousImage, alg=NFFTAlg())
    cache = create_cache(alg, parent(model), model.kernel)
    return ModelImage(model, parent(model), cache)
end

"""
    modelimage(img::ContinuousImage, cache::AbstractCache)

Create a model image directly using an image, i.e. treating it as the model. Additionally
reuse a previously compute image `cache`. This can be used when directly modeling an
image of a fixed size and number of pixels.
"""
@inline function modelimage(img::ContinuousImage, cache::AbstractCache)
    newcache = update_cache(cache, parent(img), img.kernel)
    return ModelImage(img, parent(img), newcache)
end
