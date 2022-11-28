export ContinuousImage

"""
    ContinuousImage{A<:AbstractDimArray, P}


"""
struct ContinuousImage{A <: KeyedArray, P} <: AbstractModel
    """
    Discrete representation of the image. This must be a DimArray where at least two of the
    """
    img::A
    """
    Image Kernel that transforms from the discrete image to a continuous one. This is
    sometimes called a pulse function in `eht-imaging`.
    """
    kernel::P
end

function ContinuousImage(im::AbstractMatrix, fovx::Real, fovy::Real, x0::Real, y0::Real, pulse, header=nothing)
    xitr, yitr = imagepixels(fovx, fovy, size(img, 1), size(img,2), x0, y0)
    img = IntensityMap(im, (X=xitr, Y=yitr), header)
    return ContinuousImage(img, pulse)
end

function ContinuousImage(im::AbstractMatrix, fov::Real, x0::Real, y0::Real, pulse, header=nothing)
    return ContinuousImage(im, fov, fov, x0, y0, pulse, header)
end

Base.getindex(img::ContinuousImage, args...) = getindex(img.img, args...)
Base.setindex!(img::ContinuousImage, args...) = setindex!(img.img, args...)


ComradeBase.imagepixels(img::ContinuousImage) = NamedTuple{names.(dims(img.img))}(dims(img.img))

# IntensityMap will obey the Comrade interface. This is so I can make easy models
visanalytic(::Type{<:ContinuousImage}) = NotAnalytic() # not analytic b/c we want to hook into FFT stuff
imanalytic(::Type{<:ContinuousImage}) = IsAnalytic()
isprimitive(::Type{<:ContinuousImage}) = IsPrimitive()

function intensity_point(m::ContinuousImage, p)
    dx, dy = pixelsizes(m)
    sum = zero(eltype(m.img))
    @inbounds for (I, p0) in pairs(grid(m.img))
        dp = (p.X - p0.X, p.Y - p0.Y)
        k = intensity_point(m.pulse, dp)
        sum += m.img[I]*k/(dx*dy)
    end
    return sum
end
