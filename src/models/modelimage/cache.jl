export create_cache, update_cache, DFTAlg

abstract type FourierTransform end

abstract type NUFT <: FourierTransform end

"""
    $(TYPEDEF)
This defines an abstract cache that can be used to
hold or precompute some computations.
"""
abstract type AbstractCache end


"""
    create_cache(alg, img)
Creates a Fourier transform cache for a img. This is usually called internally for ModelImage.
However, in certain circumstances the end-user can call this to speed up a method.
"""
function create_cache end

"""
    update_cache(cache, img)
Update the Fourier transform cache. This will reuse an FFT/NFFT plan saving some computational
time.
"""
function update_cache end

"""
    $(TYPEDEF)
Fourier transform type that specifies we will use
the FFTW package to compute the Fourier transform.

# Fields
$(FIELDS)

"""
Base.@kwdef struct FFTAlg <: FourierTransform
    """
    The amount to pad the image by.
    Note we actually round up to the nearest factor
    of 2, but this will be improved in the future to use
    small primes
    """
    padfac::Int = 2
end
include(joinpath(@__DIR__, "fft_alg.jl"))

struct ObservedNUFT{A<:NUFT, T} <: NUFT
    alg::A
    uv::T
end
padfac(a::ObservedNUFT) = padfac(a.alg)

struct NUFTCache{A,P,M,PI,I} <: AbstractCache
    alg::A
    plan::P
    phases::M
    pulse::PI
    img::I
end
include(joinpath(@__DIR__, "nuft.jl"))
