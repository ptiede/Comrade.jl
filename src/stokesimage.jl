export IntensityMap, PolarizedMap, fov, imagepixels, pixelsizes,
       SqExpKernel, BSplineKernel, stokes_parameter

abstract type AbstractIntensityMap{T,S} <: AbstractMatrix{T} end
abstract type AbstractPolarizedMap end

function fov(im::AbstractIntensityMap)
    return im.fovx, im.fovy
end

struct PolarizedMap{SI<:AbstractIntensityMap,
                    SQ<:AbstractIntensityMap,
                    SU<:AbstractIntensityMap,
                    SV<:AbstractIntensityMap} <: AbstractPolarizedMap
    I::SI
    Q::SQ
    U::SU
    V::SV
    function PolarizedMap(I::SI,Q::SQ,U::SU,V::SV) where {SI, SQ, SU, SV}
        @assert size(I) == size(Q) == size(U) == size(V) "Image sizes must be equal in polarized map"
        @assert fov(I) == fov(Q) == fov(U) == fov(V) "Image fov must be equal in polarized map"
        new{SI,SQ,SU,SV}(I,Q,U,V)
    end
end



Base.Base.@propagate_inbounds function Base.getindex(pimg::PolarizedMap, i...)
    return StokesVector(pimg.I[i...], pimg.Q[i...], pimg.U[i...], pimg.V[i...])
end

@inline stokes_parameter(pimg::PolarizedMap, p::Symbol) = getproperty(pimg, p)


# Base.IndexStyle(::Type{<:AbstractPolarizedMap{T,S}}) where {T, X, S, K} = Base.IndexStyle(S)
# Base.size(im::PolarizedMap) = size(im.pim)
# Base.getindex(im::PolarizedMap, i::Int) = getindex(im.pim, i)
# Base.getindex(im::PolarizedMap, I...) = getindex(im.pim, I...)
# Base.setindex!(im::PolarizedMap, x, i::Int) = setindex!(im.pim, x, i)


# function PolarizedMap(, fovx, fovy, pulse=BSplineKernel{0}())
#     @assert size(pim,3) == 4 "There must be 4 polarizations"
#     ny,nx,_ = size(pim)
#     psizex=fovx/(nx-1)
#     psizey=fovy/(ny-1)
#     return PolarizedMap(pim,
#                        convert(typeof(psizex),fovx),
#                        convert(typeof(psizey),fovy),
#                        psizex,
#                        psizey,
#                        pulse)

# end


# function Base.getproperty(pim::PolarizedMap, x::Symbol)
#     if x ∈ propertynames(pim)
#         return getfield(pim, x)
#     end
#     ind = findfirst(==(x), (:I,:Q,:U,:V))
#     isnothing(ind) && throw("type PolarizedMap has no field $(x)")
#     im = @view pim[:,:,ind]
#     return IntensityMap(im, getfield(pim, :fovx), getfield(pim, :fovy), getfield(pim, :pulse))
# end



"""
    κflux(ImageKernel)
Finds the total flux for the image kernel. This is important for find the
flux of the model.
"""
function κflux end

"""
    κ(ImageKernel, x)
Find the pointwise value of the image interpolation kernel `κ` at `x`.
"""
function κ end

"""
    ImageKernel
Pixel response function for a radio image model. This makes
a discrete sampling continuous by picking a certain *smoothing*
kernel for the image.

# Notes
To see the implemented ImageKernels please use the subtypes function i.e.
`subtypes(ImageKernel)`
"""
abstract type ImageKernel <: AbstractModel end

visanalytic(::Type{<:ImageKernel}) = IsAnalytic()
imanalytic(::Type{<:ImageKernel}) = IsAnalytic()
isprimitive(::Type{<:ImageKernel}) = IsPrimitive()

struct DeltaKernel{T} <: ImageKernel end
DeltaKernel() = DeltaKernel{Float64}()

# This should really be a delta function but whatever
intensity_point(::DeltaKernel{T}, x, y) = one(T)
visibility_point(::DeltaKernel{T}, x, y) = one(T)


"""
    $(TYPEDEF)
Normalized square exponential kernel, i.e. a Gaussian. Note the
smoothness is modfied with `ϵ` which is the inverse variance in units of
1/pixels².
"""
struct SqExpKernel{T} <: ImageKernel
    ϵ::T
end
@inline @fastmath κ(b::SqExpKernel, x) = exp(-0.5*b.ϵ^2*x^2)/sqrt(2*π/b.ϵ^2)
@inline κflux(::SqExpKernel{T}) where {T} = one(T)
@inline @fastmath ω(b::SqExpKernel, u) = exp(-2*(π*u/b.ϵ)^2)


@doc raw"""
    $(TYPEDEF)
Uses the basis spline (BSpline) kernel of order `N`. These are the kernel that come
from recursively convolving the tophat kernel
```math
    B_0(x) = \begin{cases} 1 & |x| < 1 \\ 0 & otherwise \end{cases}
```
`N` times.

## Notes

BSpline kernels have a number of nice properties:
1. Simple frequency response $\sinc(u/2)^N$
2. preserve total intensity

For `N`>1 these kernels aren't actually interpolation kernels however, this doesn't matter
for us.

Currently only the 0,1,3 order kernels are implemented.
"""
struct BSplineKernel{N} <: ImageKernel end
@inline ω(::BSplineKernel{N}, u) where {N} = sinc(π*u)^N
@inline κflux(::BSplineKernel) = 1.0

@inline κ(b::BSplineKernel{0}, x::T) where {T} = abs(x) < 0.5 ? one(T) : zero(T)

@inline function κ(b::BSplineKernel{1}, x::T) where {T}
    mag = abs(x)
    return mag < 1 ? 1-mag : zero(T)
end

@inline function κ(b::BSplineKernel{3}, x::T) where {T}
    mag = abs(x)
    if mag < 1
        return evalpoly(mag, (4, 0, -6, 3))/6
    elseif 1 ≤ mag < 2
        return evalpoly(mag, (8, -12, 6, -1))/6
    else
        return zero(T)
    end
end
