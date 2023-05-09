padfac(alg::NUFT) = alg.padfac

function padimage(::NUFT, img::IntensityMapTypes)
    #pf = padfac(alg)
    #cimg = convert(Matrix{Complex{eltype(img)}}, img.img)
    return img
    # ny,nx = size(img)
    # nnx = nextpow(2, pf*nx)
    # nny = nextpow(2, pf*ny)
    # nsx = nnx÷2-nx÷2
    # nsy = nny÷2-ny÷2
    # cimg = convert(Matrix{Complex{eltype(img)}}, img)
    # return PaddedView(zero(eltype(cimg)), cimg,
    #                   (1:nnx, 1:nny),
    #                   (nsx+1:nsx+nx, nsy+1:nsy+ny)
    #                  )
end

padimage(alg::ObservedNUFT, img::IntensityMapTypes) = padimage(alg.alg, img)


function create_cache(alg::ObservedNUFT, img::IntensityMapTypes, pulse::Pulse=DeltaPulse())
    # pimg = padimage(alg, img)

    # make nuft plan
    plan = plan_nuft(alg, img)
    # get phases and pulse functions
    phases = make_phases(alg, img, pulse)

    return create_cache(alg, plan, phases, img, pulse)
end

function create_cache(alg::NUFT, img::IntensityMapTypes, pulse::Pulse=DeltaPulse())
    # pimg = padimage(alg, img)
    return NUFTCache(alg, nothing, nothing, pulse, img)
end

function update_cache(cache::NUFTCache, img::IntensityMapTypes, pulse::Pulse=DeltaPulse())
    # # pimg = padimage(cache.alg, img)
    # cache2 = update_phases(cache, img, pulse)
    return create_cache(cache.alg, cache.plan, cache.phases, img, pulse)
    # return cache
end

function update_phases(cache::NUFTCache, img::IntensityMapTypes, pulse::Pulse)
    #if cache.pulse != pulse
    #    phases = make_phases(cache.alg, img, pulse)
    #    return @set cache.phases = phases
    #else
        return cache
    #end
end

function nocachevis(m::ModelImage{M,I,<:NUFTCache}, u, v, time, freq) where {M,I<:IntensityMap}
    alg = ObservedNUFT(m.cache.alg, vcat(u', v'))
    cache = create_cache(alg, m.image)
    m = @set m.cache = cache
    return visibilities_numeric(m, u, v, time, freq)
end

function checkuv(uv, u, v)
    @assert u == @view(uv[1,:]) "Specified u don't match uv in cache. Did you pass the correct u,v?"
    @assert v == @view(uv[2,:]) "Specified v don't match uv in cache. Did you pass the correct u,v?"
end



#using ReverseDiff
#using NFFT
#ReverseDiff.@grad_from_chainrules nuft(A, b::ReverseDiff.TrackedArray)
#ReverseDiff.@grad_from_chainrules nuft(A, b::Vector{<:ReverseDiff.TrackedReal})


ChainRulesCore.@non_differentiable checkuv(alg, u::AbstractArray, v::AbstractArray)

function visibilities_numeric(m::ModelImage{M,I,<:NUFTCache{A}},
                      u, v, time, freq) where {M,I<:IntensityMap,A<:ObservedNUFT}
    checkuv(m.cache.alg.uv, u, v)
    vis =  nuft(m.cache.plan, complex.(Comrade.baseimage(m.cache.img)))
    return conj.(vis).*m.cache.phases
end

function visibilities_numeric(m::ModelImage{M,I,<:NUFTCache{A}},
                      u, v, time, freq) where {M,I<:StokesIntensityMap,A<:ObservedNUFT}
    checkuv(m.cache.alg.uv, u, v)
    visI =  conj.(nuft(m.cache.plan, complex.(Comrade.baseimage(stokes(m.cache.img, :I))))).*m.cache.phases
    visQ =  conj.(nuft(m.cache.plan, complex.(Comrade.baseimage(stokes(m.cache.img, :Q))))).*m.cache.phases
    visU =  conj.(nuft(m.cache.plan, complex.(Comrade.baseimage(stokes(m.cache.img, :U))))).*m.cache.phases
    visV =  conj.(nuft(m.cache.plan, complex.(Comrade.baseimage(stokes(m.cache.img, :V))))).*m.cache.phases
    r = StructArray{StokesParams{eltype(visI)}}((I=visI, Q=visQ, U=visU, V=visV))
    return r
end





function visibilities_numeric(m::ModelImage{M,I,<:NUFTCache{A}},
                      u, v, time, freq) where {M,I,A<:NUFT}
    return nocachevis(m, u, v, time, freq)
end


"""
    NFFTAlg
Uses a non-uniform FFT to compute the visibilities.
You can optionally pass uv which are the uv positions you will
compute the NFFT at. This can allow for the NFFT plan to be cached improving
performance

# Fields
$(FIELDS)

"""
Base.@kwdef struct NFFTAlg{T,N,F} <: NUFT
    """
    Amount to pad the image
    """
    padfac::Int = 1
    """
    Kernel size parameters. This controls the accuracy of NFFT you do not usually need to change this
    """
    m::Int = 4
    """
    Over sampling factor. This controls the accuracy of NFFT you do not usually need to change this.
    """
    σ::T = 2.0
    """
    Window function for the NFFT. You do not usually need to change this
    """
    window::Symbol = :kaiser_bessel
    """
    NFFT interpolation algorithm. TENSOR is the fastest but takes the longest to precompute
    """
    precompute::N=NFFT.TENSOR
    """
    Flag block partioning should be used to speed up computation
    """
    blocking::Bool = true
    """
    Flag if the node should be sorted in a lexicographic way
    """
    sortNodes::Bool = false
    """
    Flag if the deconvolve indices should be stored, Currently required for GPU
    """
    storeDeconvolutionIdx::Bool = true
    """
    Flag passed to inner AbstractFFT. The fastest FFTW is FFTW.MEASURE but takes the longest
    to precompute
    """
    fftflags::F = FFTW.MEASURE
end
include(joinpath(@__DIR__, "nfft_alg.jl"))

"""
    DFTAlg
Uses a discrete fourier transform. This is not very efficient for larger images. In those cases
 NFFTAlg or FFTAlg are more reasonable. For small images this is a reasonable choice especially
since it's easy to define derivatives.
"""
struct DFTAlg <: NUFT end
include(joinpath(@__DIR__, "dft_alg.jl"))
