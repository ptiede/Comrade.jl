"""
    $(TYPEDEF)

Abstract type that specified which fourier transform to use
"""
abstract type FourierTransform end

struct DFT <: FourierTransform end

"""
    $(TYPEDEF)
This defines an abstract cache that can be used to
hold or precompute some computations.
"""
abstract type AbstractCache end

"""
    $(TYPEDEF)
No cache is used. This is typically used when the model is
analytic in the Fourier domain.
"""
struct NoCache <: AbstractCache end


"""
    $(TYPEDEF)
Fourier transform type that specifies we will use
the FFTW package to compute the Fourier transform.

# Fields
$(FIELDS)

"""
Base.@kwdef struct FFT <: FourierTransform
    """
    The amount to pad the image by.
    Note we actually round up to the nearest factor
    of 2, but this will be improved in the future to use
    small primes
    """
    padfac::Int = 2
end

function create_interpolator(u, v, vis)
    # Construct the interpolator
    itp = interpolate(vis, BSpline(Cubic(Line(OnGrid()))))
    etp = extrapolate(itp, zero(eltype(vis)))
    scale(etp, u, v)
end

"""
    $(TYPEDEF)
The cache used when the `FFT` algorithm is used to compute
visibilties. This is an internal type and is not part of the public API

# Fields
$(FIELDS)
"""
struct FFTCache{P,I} <: AbstractCache
    """ FFTW Plan"""
    plan::P
    """FFT interpolator function"""
    sitp::I
end

"""
    $(SIGNATURES)
Creates the model cache given for the algorithm `alg`
using the `model` and a image cache `image`
"""
function create_cache(alg::FFT, model, img)
    intensitymap!(img, model)
    dx,dy = pixelsizes(img)
    ny,nx = size(img)
    padfac = alg.padfac
    nnx = nextpow(2, padfac*nx)
    nny = nextpow(2, padfac*ny)
    pimg = PaddedView(zero(eltype(img)), img.im, (nnx, nny))

    uu = fftshift(fftfreq(nnx, 1/dx))
    vv = fftshift(fftfreq(nny, 1/dy))
    #plan = plan_fft(pimg)
    vis = fftshift(fft(pimg))
    #println(sum(img)*dx*dy)
    #println(sum(pimg)*dx*dy)
    x0,y0 = first.(imagepixels(img))
    @inbounds @fastmath for I in CartesianIndices(vis)
        iy, ix = Tuple(I)
        vis[I] = conj(vis[I])*dx*dy*exp(2im*π*(uu[ix]*x0 + vv[iy]*y0))
    end
    sitp = create_interpolator(uu, vv, vis)
    return FFTCache(nothing, sitp)
end


function fouriermap(m, img)
    ny,nx = size(img)
    x,y = imagepixels(img)
    dx,dy = pixelsizes(img)

    uu = fftshift(fftfreq(length(x), 1/dx))
    vv = fftshift(fftfreq(length(y), 1/dy))

    vis = Matrix{Complex{eltype(img)}}(undef, ny, nx)

    @inbounds for I in CartesianIndices(vis)
        iy, ix = Tuple(I)
        vp = visibility(m, uu[ix], vv[iy])
        vis[I] = vp
    end
    return vis
end

function phasedecenter!(vis, img)
    ny,nx = size(img)
    x,y = imagepixels(img)
    dx,dy = pixelsizes(img)

    uu = fftshift(fftfreq(length(x), 1/dx))
    vv = fftshift(fftfreq(length(y), 1/dy))

    for I in CartesianIndices(vis)
        iy, ix = Tuple(I)
        vis[I] = conj(vis[I])*exp(-2im*π*(uu[ix]*first(x) + vv[iy]*first(y)))*nx*ny/(dx*dy)
    end
    return vis
end
