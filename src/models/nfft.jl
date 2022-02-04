function padimage(img, alg::NFFT)
    padfac = alg.padfac
    ny,nx = size(img)
    nnx = nextpow(2, padfac*nx)
    nny = nextpow(2, padfac*ny)
    nsx = nnx÷2-nx÷2
    nsy = nny÷2-ny÷2
    return PaddedView(zero(eltype(img)), img.im,
                      (1:nnx, 1:nny),
                      (nsx+1:nsx+nx, nsy+1:nsy+ny)
                     )
end

function create_cache(alg::NFFT, img)
    pimg = padimage(img, alg)

    # No supplied uv positions so can't plan NFFT
    isnothing(alg.uv) && return NFFTCache(alg, nothing, LinearAlgebra.I, pimg)

    # User supplied
    plan = plan_nfft(alg.uv, size(pimg), alg.m)
    phases = @. exp(1im*π*(u*dx + v*dy))*dx*dy#*exp(1im*π*(u*x + v*y))*dx*dy
    return plan, phases, pimg
end

Base.@kwdef struct NFFT <: FourierTransform
    padfac::Int = 2
    m::Int = 10
    uv::T = nothing
end

struct NFFTCache{A<:NFFT,C,M,I} <: AbstractCache
    alg::A
    plan::P
    phases::M
    pimg::I
end
