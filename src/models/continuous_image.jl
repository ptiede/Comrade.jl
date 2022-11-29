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

function ContinuousImage(img::IntensityMap, pulse::AbstractModel)
    dx, dy = pixelsizes(img)
    spulse = stretched(pulse, dx, dy)
    return ContinuousImage{typeof(img), typeof(spulse)}(img, spulse)
end


function ContinuousImage(im::AbstractMatrix, fovx::Real, fovy::Real, x0::Real, y0::Real, pulse, header=nothing)
    xitr, yitr = imagepixels(fovx, fovy, size(img, 1), size(img,2), x0, y0)
    img = IntensityMap(im, (X=xitr, Y=yitr), header)
    spulse = stretched(pulse, step(xitr), step(yitr))
    return ContinuousImage(img, spulse)
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

radialextent(c::ContinuousImage) = maximum(values(fieldofview(c.img)))/2

function intensity_point(m::ContinuousImage, p)
    sum = zero(eltype(m.img))
    @inbounds for (I, p0) in pairs(grid(m.img))
        dp = (X=(p.X - p0.X), Y=(p.Y - p0.Y))
        k = intensity_point(m.kernel, dp)
        sum += m.img[I]*k
    end
    return sum
end

convolved(cimg::ContinuousImage, m::AbstractModel) = ContinuousImage(cimg.img, convolved(cimg.pulse, m))
