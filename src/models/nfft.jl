export NFFTAlg

Base.@kwdef struct NFFTAlg{T} <: FourierTransform
    padfac::Int = 1
    m::Int = 10
    uv::T = nothing
end

function NFFTAlg(u::AbstractArray, v::AbstractArray; padfac=2, m=10)
    uv = Matrix{eltype(u)}(undef, 2, length(u))
    uv[1,:] .= u
    uv[2,:] .= v
    return NFFTAlg(;padfac, m, uv)
end

function NFFTAlg(ac::ArrayConfiguration; padfac=2, m=10)
    u, v = getuv(ac)
    return NFFTAlg(u, v, padfac, m=10)
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

@fastmath function create_cache(alg::NFFTAlg, img)
    pimg = padimage(img, alg)

    # No supplied uv positions so can't plan NFFT
    isnothing(alg.uv) && return NFFTCache(alg, nothing, LinearAlgebra.I, pimg)
    dx, dy = pixelsizes(img)
    # User supplied
    uv2 = copy(alg.uv)
    uv2[1,:] .= uv2[1,:]*dx
    uv2[2,:] .= uv2[2,:]*dx
    plan = plan_nfft(uv2, size(pimg); precompute=NFFT.POLYNOMIAL)
    phases = @. cispi((@view(uv2[1,:]) + @view(uv2[2,:])))*dx*dy
    return NFFTCache(alg, plan, phases, pimg)
end

function phasenfft!(vis, phases)
    for I in eachindex(vis, phases)
        vis[I] = conj(vis[I])*phases[I]
    end
    return vis
end

function visibilities(cache::NFFTCache)
    vis = cache.plan*cache.pimg'
    return phasenfft!(vis, cache.phases)
end

function nocachevis(m::ModelImage, u, v)
    alg = NFFTAlg(u,v; padfac=m.cache.alg.padfac, m=m.cache.alg.m)
    cache = create_cache(alg, m.image)
    m = @set m.cache = cache
    visibilities(m, u, v)
end

function checkuv(uv, u, v)
    @assert u == @view(uv[1,:]) "Specified u don't match uv in cache. Did you pass the correct u,v to NFFTAlg?"
    @assert v == @view(uv[2,:]) "Specified u don't match uv in cache. Did you pass the correct u,v to NFFTAlg?"
end


function visibilities(m::ModelImage{M,I,C}, u::AbstractArray, v::AbstractArray) where {M,I,C<:NFFTCache}
    if isnothing(m.cache.plan)
        return nocachevis(m, u, v)
    else
        checkuv(m.cache.alg.uv, u, v)
        vis =  m.cache.plan*m.cache.pimg
        vis .= vis.*m.cache.phases
        return vis
    end
end

function amplitudes(m::ModelImage{M,I,C}, u::AbstractArray, v::AbstractArray) where {M,I,C<:NFFTCache}
    if isnothing(m.cache.plan)
        return nocachevis(m, u, v)
    else
        checkuv(m.cache.alg.uv, u, v)
        vis =  m.cache.plan*m.cache.pimg
        vis .= abs.(vis)
        return vis
    end
end

function _recache(m)
    c = NFFTCache(m.cache.alg, nothing, nothing, m.cache.pimg)
    m2 = ModelImage(m.model, m.image, c)
    return m2
end

function bispectra(m::ModelImage{M,I,C},
                    u1::AbstractArray,
                    v1::AbstractArray,
                    u2::AbstractArray,
                    v2::AbstractArray,
                    u3::AbstractArray,
                    v3::AbstractArray) where {M,I,C<:NFFTCache}
    # TODO: Fix this so we can actually use a stored cache for closures
    n = length(u1)
    # In testing it is faster to concatenate everything into 1 vector and then
    # take an NFFT of that
    u = vcat(u1,u2,u3)
    v = vcat(v1,v2,v3)
    vis = nocachevis(m, u, v)
    @inbounds (@view(vis[1:n])).*(@view(vis[n+1:2n])).*(@view vis[2n+1:3n])
end

function closure_phases(m::ModelImage{M,I,C},
                        u1::AbstractArray,
                        v1::AbstractArray,
                        u2::AbstractArray,
                        v2::AbstractArray,
                        u3::AbstractArray,
                        v3::AbstractArray) where {M,I,C<:NFFTCache}
    # TODO: Fix this so we can actually use a stored cache for closures
    n = length(u1)
    # In testing it is faster to concatenate everything into 1 vector and then
    # take an NFFT of that
    u = vcat(u1,u2,u3)
    v = vcat(v1,v2,v3)
    vis = nocachevis(m, u, v)
    phase = @inbounds angle.(@view(vis[1:n])) .+ angle.(@view(vis[n+1:2n])) .+ angle.(@view vis[2n+1:3n])
    return phase
end

function logclosure_amplitudes(m::ModelImage{M,I,C},
                               u1::AbstractArray,
                               v1::AbstractArray,
                               u2::AbstractArray,
                               v2::AbstractArray,
                               u3::AbstractArray,
                               v3::AbstractArray,
                               u4::AbstractArray,
                               v4::AbstractArray,
                               ) where {M,I,C<:NFFTCache}
    # TODO: Fix this so we can actually use a stored cache for closures
    n = length(u1)
    u = vcat(u1,u2,u3,u4)
    v = vcat(v1,v2,v3,v4)
    vis = nocachevis(m, u, v)
    amp1 = @view vis[1:n]
    amp2 = @view vis[n+1:2n]
    amp3 = @view vis[2n+1:3n]
    amp4 = @view vis[3n+1:4n]

    lcamp = @. log(amp1*amp2/(amp3*amp4))
    return lcamp
end
