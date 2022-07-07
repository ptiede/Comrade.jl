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
    create_cache(alg::AbstractFourierTransform, img::AbstractIntensityMap)

Creates a Fourier transform cache for a img using algorithm `alg`. For non-analytic visibility models this
can significantly speed up computations.

# Examples

```julia-repl
julia> u,v = rand(100), rand(100)
julia> cache = create_cache(DFTAlg(u, v), IntensityMap(randn(50,50), 10.0, 10.0))
julia> cache = create_cache(NFFTAlg(u, v), IntensityMap(randn(50,50), 10.0, 10.0))
julia> cache = create_cache(FFTAlg(), IntensityMap(randn(50,50), 10.0, 10.0))
```
"""
function create_cache end

"""
    update_cache(cache, img)

Update the Fourier transform cache. This will reuse an FFT/NFFT plan saving some computational
time.

# Note

This is an intenal method than an end user shouldn't have to usually call.
"""
function update_cache end

"""
    $(TYPEDEF)

Use an FFT to compute the approximate numerical visibilities of a model.
For a DTFT see [`DFTAlg`](@ref DFTAlg) or for an NFFT [`NFFTAlg`](@ref NFFTAlg)

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

"""
    $(TYPEDEF)

Container type for a non-uniform Fourier transform (NUFT).
This stores the uv-positions that the model will be sampled at in the Fourier domain,
allowing certain transformtion factors (e.g., NUFT matrix) to be cached.

This is an internal type, an end user should instead create this using [`NFFTAlg`](@ref NFFTAlg)
or [`DFTAlg`](@ref DFTAlg).
"""
struct ObservedNUFT{A<:NUFT, T} <: NUFT
    """
    Which NUFT algorithm to use (e.g. NFFTAlg or DFTAlg)
    """
    alg::A
    """
    uv positions of the NUFT transform. This is used for precomputation.
    """
    uv::T
end
padfac(a::ObservedNUFT) = padfac(a.alg)

"""
    $(TYPEDEF)

Internal type used to store the cache for a non-uniform Fourier transform (NUFT).

The user should instead create this using the [`create_cache`](@ref create_cache) function.
"""
struct NUFTCache{A,P,M,PI,I} <: AbstractCache
    alg::A # which algorithm to use
    plan::P #NUFT matrix or plan
    phases::M #FT phases needed to phase center things
    pulse::PI #pulse function to make continuous image
    img::I # image buffer to hold the image of the model that we will take a FT of.
end
include(joinpath(@__DIR__, "nuft.jl"))
