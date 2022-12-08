padfac(alg::NUFT) = alg.padfac

function padimage(::NUFT, img::Union{StokesIntensityMap, IntensityMap})
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

padimage(alg::ObservedNUFT, img::Union{StokesIntensityMap,IntensityMap}) = padimage(alg.alg, img)



function create_cache(alg::ObservedNUFT, img::Union{StokesIntensityMap,IntensityMap}, pulse=DeltaPulse())
    pimg = padimage(alg, img)

    # make nuft plan
    dx, dy = pixelsizes(img)
    plan = plan_nuft(alg, pimg)
    # get phases and pulse functions
    phases = make_phases(alg, img, pulse)

    return create_cache(alg, plan, phases, pimg, pulse)
end

function create_cache(alg::NUFT, img::Union{StokesIntensityMap,IntensityMap}, pulse=DeltaPulse())
    pimg = padimage(alg, img)
    return NUFTCache(alg, nothing, nothing, pimg, pulse)
end

function update_cache(cache::NUFTCache, img::Union{StokesIntensityMap,IntensityMap}, pulse=DeltaPulse())
    pimg = padimage(cache.alg, img)
    cache2 = update_phases(cache, img, pulse)
    create_cache(cache2.alg, cache2.plan, cache2.phases, pimg, pulse)
end

function update_phases(cache::NUFTCache, img::IntensityMap, pulse)
    #if cache.pulse != pulse
    #    phases = make_phases(cache.alg, img, pulse)
    #    return @set cache.phases = phases
    #else
        return cache
    #end
end

function nocachevis(m::ModelImage{M,I,<:NUFTCache}, p) where {M,I<:IntensityMap}
    u = p.U
    v = p.V
    alg = ObservedNUFT(m.cache.alg, vcat(u', v'))
    cache = create_cache(alg, m.image)
    m = @set m.cache = cache
    return _visibilities(m, StructArray((U=u, V=v)))
end

function checkuv(uv, u, v)
    @assert u == @view(uv[1,:]) "Specified u don't match uv in cache. Did you pass the correct u,v?"
    @assert v == @view(uv[2,:]) "Specified v don't match uv in cache. Did you pass the correct u,v?"
end


function nuft(A, b)
    return A*b
end

function ChainRulesCore.rrule(::typeof(nuft), A::NFFTPlan, b)
    #pr = ChainRulesCore.ProjectTo(b)
    vis = A*b
    function nuft_pullback(Δy)
        Δf = NoTangent()
        dy = similar(vis)
        dy .= unthunk(Δy)
        ΔA = A'*dy
        return Δf, NoTangent(), ΔA
    end
    return vis, nuft_pullback
end

#using ReverseDiff
#using NFFT
#ReverseDiff.@grad_from_chainrules nuft(A, b::ReverseDiff.TrackedArray)
#ReverseDiff.@grad_from_chainrules nuft(A, b::Vector{<:ReverseDiff.TrackedReal})


ChainRulesCore.@non_differentiable checkuv(alg, u::AbstractArray, v::AbstractArray)

function _visibilities(m::ModelImage{M,I,<:NUFTCache{A}},
                      p) where {M,I<:IntensityMap,A<:ObservedNUFT}
    u = p.U
    v = p.V
    checkuv(m.cache.alg.uv, u, v)
    vis =  nuft(m.cache.plan, complex.(m.cache.img))
    return conj.(vis).*m.cache.phases
end

function _visibilities(m::ModelImage{M,I,<:NUFTCache{A}},
                      p) where {M,I<:StokesIntensityMap,A<:ObservedNUFT}
    u = p.U
    v = p.V
    checkuv(m.cache.alg.uv, u, v)
    visI =  nuft(m.cache.plan, complex.(stokes(m.cache.img, :I)))
    visQ =  nuft(m.cache.plan, complex.(stokes(m.cache.img, :Q)))
    visU =  nuft(m.cache.plan, complex.(stokes(m.cache.img, :U)))
    visV =  nuft(m.cache.plan, complex.(stokes(m.cache.img, :V)))
    return StructArray{StokesParams{eltype(visI)}}((I=visI, Q=visQ, U=visU, V=visV)).*m.cache.phases
end



function _visibilities(m::ModelImage{M,I,<:NUFTCache{A}},
                      p) where {M,I,A<:NUFT}
    return nocachevis(m, p)
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
    Flag blcok partioning should be used to speed up computation
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
