export centroid_mean, center_image, convolve!, convolve

function centroid_mean(imgs::AbstractVector{<:IntensityMap})
    return mapreduce(+, imgs) do img
        center_image(img)
    end
    return mimg./length(imgs)
end

"""
    center_image(img::IntensityMap)

centers the `img` such that the centroid of the image is approximately at the origin.
"""
function center_image(img::IntensityMap)
    x, y = centroid(img)
    return modify(img, Shift(-x, -y))
end

# This is an internal struct that is use to modify IntensityMaps so that we can hook into
# Comrade image modifier interface.
struct InterpolatedImage{I,P} <: AbstractModel
    img::I
    itp::P
    function InterpolatedImage(img::IntensityMap)
        itp = BicubicInterpolator(img.X, img.Y, img, StrictBoundaries())
        return new{typeof(img), typeof(itp)}(img, itp)
    end
end


imanalytic(::Type{<:InterpolatedImage})  = IsAnalytic()
visanalytic(::Type{<:InterpolatedImage}) = NotAnalytic()

function intensity_point(m::InterpolatedImage, p)
    (m.img.X[begin] > p.X || p.X > m.img.X[end]) && return zero(p.X)
    (m.img.Y[begin] > p.Y || p.Y > m.img.Y[end]) && return zero(p.X)
    return m.itp(p.X, p.Y)/(step(m.img.X)*step(m.img.Y))
end
function ModifiedModel(img::IntensityMap, transforms)
    ms = ModifiedModel(InterpolatedImage(img), transforms)
    return intensitymap(ms, axiskeys(img))
end


"""
    modify(img::IntensityMap, transforms...)

This modifies the `img` by applying the `transforms...` returning a transformed `IntensityMap`

!!! note
Unlike when `modify` is applied to a `<:AbstractModel` this returns an already modified image.
"""
modify(img::IntensityMap, transforms...) = ModifiedModel(img, transforms)


"""
    convolve!(img::IntensityMap, m::AbstractModel)

Convolves an `img` with a given analytic model `m`. This is useful for blurring the
image with some model. For instance to convolve a image with a Gaussian you would do
```julia
convolve!(img, Gaussian())
```

# Notes
This method does not automatically pad your image. If there is substantial flux at the boundaries
you will start to see artifacts.
"""
function convolve!(img::IntensityMap{<:Real}, m::AbstractModel)
    @assert visanalytic(typeof(m)) isa IsAnalytic "Convolving model must have an analytic Fourier transform currently"
    p = plan_rfft(baseimage(img))

    (;X, Y) = imagepixels(img)
    # plan_rfft uses just the positive first axis to respect real conjugate symmetry
    U = rfftfreq(size(img, 1), inv(step(X)))
    V = fftfreq(size(img, 2), inv(step(Y)))

    # TODO maybe ask a user to pass a vis buffer as well?
    vis = p*baseimage(img)

    # Conjugate because Comrade uses +2πi exponent
    vis .*= conj(visibility_point.(Ref(m), U, V', 0, 0))
    pinv = plan_irfft(vis, size(img, 1))
    mul!(baseimage(img), pinv, vis)
    return img
end

"""
    convolve(img::IntensityMap, m::AbstractModel)

Convolves an `img` with a given analytic model `m`. This is useful for blurring the
image with some model. For instance to convolve a image with a Gaussian you would do
```julia
convolve(img, Gaussian())
```

For the inplace version of the function see [`convolve!`](@ref)

# Notes
This method does not automatically pad your image. If there is substantial flux at the boundaries
you will start to see artifacts.
"""
function convolve(img::IntensityMap, m::AbstractModel)
    cimg = copy(img)
    return convolve!(cimg, m)
end

function convolve!(img::IntensityMap{<:StokesParams}, m)
    convolve!(stokes(img, :I), m)
    convolve!(stokes(img, :Q), m)
    convolve!(stokes(img, :U), m)
    convolve!(stokes(img, :V), m)
    return img
end
