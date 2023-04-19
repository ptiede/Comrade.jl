export centroid_mean, center_image

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

"""
    modify(img::IntensityMap, transforms...)

This modifies the `img` by applying the `transforms...` returning a transformed `IntensityMap`

!!! note
Unlike when `modify` is applied to a `<:AbstractModel` this returns a already modified image.
"""
function modify(img::IntensityMap, transforms...)
    ms = modify(InterpolatedImage(img), transforms...)
    return intensitymap(ms, axiskeys(img))
end
