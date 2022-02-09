export NFFTAlg

Base.@kwdef struct NFFTAlg{T} <: FourierTransform
    padfac::Int = 2
    m::Int = 10
    uv::T = nothing
end

function NFFTAlg(u::AbstractArray, v::AbstractArray; padfac=2, m=10)
    uv = Matrix{eltype(u)}(undef, 2, length(u))
    uv[1,:] .= u
    uv[2,:] .= v
    return NFFTAlg(;padfac, m, uv)
end

struct NFFTCache{A<:NFFTAlg,P,M,I} <: AbstractCache
    alg::A
    plan::P
    phases::M
    pimg::I
end

function padimage(img, alg::NFFTAlg)
    padfac = alg.padfac
    ny,nx = size(img)
    nnx = nextpow(2, padfac*nx)
    nny = nextpow(2, padfac*ny)
    nsx = nnx÷2-nx÷2
    nsy = nny÷2-ny÷2
    cimg = convert(Matrix{Complex{eltype(img)}}, img)
    return PaddedView(zero(eltype(cimg)), cimg,
                      (1:nnx, 1:nny),
                      (nsx+1:nsx+nx, nsy+1:nsy+ny)
                     )
end

function create_cache(alg::NFFTAlg, img)
    pimg = padimage(img, alg)

    # No supplied uv positions so can't plan NFFT
    isnothing(alg.uv) && return NFFTCache(alg, nothing, LinearAlgebra.I, pimg)
    dx, dy = pixelsizes(img)
    # User supplied
    plan = plan_nfft(uv2, size(pimg))
    phases = @. cispi((@view(uv[1,:])*dx + @view(uv[2,:])*dy))*dx*dy#*exp(1im*π*(u*x + v*y))*dx*dy
    return NFFTCache(alg, plan, phases, pimg)
end

function phasenfft!(vis, phases)
    for I in eachindex(vis, phases)
        vis[I] = conj(vis[I])*phases[I]
    end
    return vis
end

function nocachevis(m::ModelImage, u, v)
    img = m.cache.pimg

    dx,dy = pixelsizes(m.image)
    uv = zeros(2, length(u))
    uv[1,:] = u*dx
    uv[2,:] = v*dy
    #compute the nfft
    #convert to complex
    p = plan_nfft(uv, size(img))
    vis = p*img
    # phase center
    phases = @. cispi((@view(uv[1,:]) + @view(uv[2,:])))*dx*dy
    phasenfft!(vis, phases)
end

function checkuv(uv, u, v)
    @assert u == @view(uv[1,:]) "Specified u don't match uv in cache. Did you pass the correct u,v to NFFTAlg?"
    @assert v == @view(uv[2,:]) "Specified u don't match uv in cache. Did you pass the correct u,v to NFFTAlg?"
end


function _visibilities(m::ModelImage{M,I,C}, u::AbstractArray, v::AbstractArray) where {M,I,C<:NFFTCache}
    if isnothing(m.cache.plan)
        return nocachevis(m, u, v)
    else
        checkuv(m.cache.alg.uv, u, v)
        vis =  m.cache.plan*m.cache.pimg
        vis .= vis.*m.cache.phases
        return vis
    end
end
