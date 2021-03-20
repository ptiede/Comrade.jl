import Base.Math.@horner

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
abstract type ImageKernel end


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
        return @horner(mag, 4, 0, -6, 3)/6
    elseif 1 ≤ mag < 2
        return @horner(mag, 8, -12, 6, -1)/6
    else
        return zero(T)
    end
end

#struct CubicKernel{T} <: ImageKernel
#    b::T
#end
#CubicKernel() = CubicKernel(-0.5)

#@inline κflux(::CubicSplineKernel) = 1.0
#@inline function κ(b::CubicSplineKernel, x::T) where {T}
#    mag = abs(x)
#    if mag ≤ 1
#        return @horner(mag, one(T), zero(T), -(b.b+T(3)), (b.b+T(2)))
#    elseif 1 < mag ≤ 2
#return b.b*@horner(mag, T(-4.0), T(8.0), T(-5.0), T(1.0))
#else
#        zero(T)
#   end
#end

#@inline @fastmath function ω(b::CubicSplineKernel, u::T) where {T}#

#end



@doc raw"""
    $(TYPEDEF)
An image model given by a set of coefficients and a kernel response or basis function.
This corresponds to a continous image defined by a finite set of points. The defined
intensity is given by
```math
    I(x,y) = \sum_{ij} c_{ij}κ(x-x_i)κ(y-y_i).
```
An important thing to note is that the ``c_{ij}`` do not represent pixel intensities, i.e.
the κ doesn't have to be an interpolating kernel.

## Example
```julia
samples = rand(10,10)
model = RImage(samples, BSplineKernel{3})
```

## Notes

This is defined in terms of pixel response, so the image size is 1μas. To resize the image
use the scale function like with other models.

## Fields

$(FIELDS)

"""
struct RImage{S,B<:ImageKernel,M<:AbstractMatrix{S}} <: AbstractRadioImage{S}
    """ Image coefficients cᵢⱼ in expansion """
    coeff::M
    """ Image kernel/basis κ that defined the delta image response """
    kernel::B
    # pixel size in 1/pixels
    psizex::S
    # pixel size in 1/pixels
    psizey::S
    function RImage(coeff::M, basis::B) where {S,M<:AbstractMatrix{S},B}
        ny, nx = size(coeff)
        psizex = one(S)/max(nx-1, 1)
        psizey = one(S)/max(ny-1, 1)
        new{S,B,M}(coeff, basis, psizex, psizey)
    end
end
FourierStyle(::Type{RImage{S,B,M}}) where {S,B,M} = IsAnalytic()

struct FourierCache{C} <: ObservationCache
    cache::C
end

function FourierCache(rimage::I, obs::Observation) where {S,I<:AbstractRadioImage{S}}
    cache = zeros(Complex{S}, size(rimage)..., nsamples(obs))
    ny,nx = size(rimage)
    u = getdata(obs, :u)
    v = getdata(obs, :v)
    dx = 1/max(nx-1,1)
    dy = 1/max(ny-1,1)
    startx = -0.5
    starty = -0.5
    x = range(startx, length=nx, step=dx)
    y = range(starty, length=ny, step=dy)
    for i in eachindex(u,v)
        cache[:,:,i] .= exp.(2im*π*(u[i].*x' .+ v[i].*y))
    end
    FourierCache(cache)
end



@inline function flux(model::RImage{S,B,M}) where {S,B,M}
    sum = zero(S)
    @avx for i in eachindex(model.coeff)
        sum += model.coeff[i]
    end
    # Divide by pixel number to convert properly
    return sum*κflux(model.kernel)/prod(size(model))
end




"""
    $(SIGNATURES)
return the size of the coefficient matrix for `model`.
"""
@inline Base.size(model::RImage) = size(model.coeff)

@inline function intensity(model::RImage{S,M,B}, x, y, args...) where {S,M,B}
    sum = zero(S)
    ny,nx = size(model)
    dx = 1/(max(nx-1,1))
    dy = 1/(max(ny-1,1))
    #The kernel is written in terms of pixel number so we convert x to it
    @inbounds @simd for I in CartesianIndices(model.coeff)
        iy,ix = Tuple(I)
        xx = x - (-0.5 + dx*(ix-1))
        yy = y - (-0.5 + dy*(iy-1))
        sum += model.coeff[I]* κ(model.kernel, xx/dx)*κ(model.kernel, yy/dy)
    end
    # Note this will be intensity per uas
    return sum
end

function cache(model::AbstractModifier, u, v)
    return cache(basemodel(model), u, v)
end

function cache(model::RImage{S,M,B}, u, v) where {S,M,B}
    ny,nx = size(model)
    dx = 1/max(nx-1,1)
    dy = 1/max(ny-1,1)
    startx = -0.5
    starty = -0.5
    upx = u*dx
    vpx = v*dy
    phasecenter = exp(2im*π*(u*startx + v*starty))
    c = zeros(Complex{S}, size(model.coeff))
    @inbounds for i in axes(model.coeff,2), j in axes(model.coeff,1)
        c[j,i] = exp(2im*π*(upx*(i-1) + vpx*(j-1)))*phasecenter*dx*dy
    end
    return c
end

@inline function visibility(model::RImage{S,M,B}, u, v, cache, args...) where {S,M,B}
    sum = zero(Complex{S})
    ny,nx = size(model)
    dx = 1/max(nx-1,1)
    dy = 1/max(ny-1,1)
    #startx = -0.5
    #starty = -0.5
    #upx = u*dx
    #vpx = v*dy
    #phasecenter = exp(2im*π*(u*startx + v*starty))
    @avx for i in axes(model.coeff,2), j in axes(model.coeff,1)
            sum += model.coeff[j,i]*cache[j,i]
    end
    return sum*dx*dy*ω(model.kernel, u*dx)*ω(model.kernel, v*dy)#*phasecenter
end
