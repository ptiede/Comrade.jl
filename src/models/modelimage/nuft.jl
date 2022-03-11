padfac(alg::NUFT) = alg.padfac

function padimage(::NUFT, img)
    #pf = padfac(alg)
    return convert(Matrix{Complex{eltype(img)}}, img)
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

padimage(alg::ObservedNUFT, img) = padimage(alg.alg, img)



function create_cache(alg::ObservedNUFT, img)
    pimg = padimage(alg.alg, img)
    # make nuft plan
    dx, dy = pixelsizes(img)
    plan = plan_nuft(alg, pimg, dx, dy)
    # get phases and pulse functions
    phases = make_phases(alg, img)

    return create_cache(alg, plan, phases, pimg)
end

function create_cache(alg::NUFT, img)
    pimg = padimage(alg, img)
    return NUFTCache(alg, nothing, nothing, pimg)
end

function update_cache(cache::NUFTCache, img)
    pimg = padimage(cache.alg, img)
    create_cache(cache.alg, cache.plan, cache.phases, pimg)
end

function nocachevis(m::ModelImage{M,I,<:NUFTCache}, u::AbstractArray, v::AbstractArray) where {M,I}
    alg = ObservedNUFT(m.cache.alg, vcat(u', v'))
    cache = create_cache(alg, m.image)
    m = @set m.cache = cache
    return visibilities(m, u, v)
end

function checkuv(uv, u, v)
    @assert u == @view(uv[1,:]) "Specified u don't match uv in cache. Did you pass the correct u,v?"
    @assert v == @view(uv[2,:]) "Specified v don't match uv in cache. Did you pass the correct u,v?"
end
using ReverseDiff


function nuft(A, b)
    return A*b
end

# function ChainRulesCore.rrule(::typeof(nuft), A, b)
#     vis = A*b
#     println("Pld")
#     function nuft_pullback(Δy)
#         Δf = NoTangent()
#         ΔA = @thunk(A'*Δy)
#         Δb = @thunk(Δy*b')
#         return Δf, ΔA, Δb
#     end
#     return vis, nuft_pullback
# end


# ReverseDiff.@grad_from_chainrules nuft(A::ReverseDiff.TrackedArray, b::ReverseDiff.TrackedArray)
# ReverseDiff.@grad_from_chainrules nuft(A, b::ReverseDiff.TrackedArray)
# ReverseDiff.@grad_from_chainrules nuft(A::ReverseDiff.TrackedArray, b)
# ReverseDiff.@grad_from_chainrules nuft(A, b::Vector{<:ReverseDiff.TrackedReal})



function _visibilities(m::ModelImage{M,I,<:NUFTCache{A}},
                      u::AbstractArray,
                      v::AbstractArray) where {M,I,A<:ObservedNUFT}
    checkuv(m.cache.alg.uv, u, v)
    #vr = real.(m.cache.plan)
    #vi = imag.(m.cache.plan)
    vis =  nuft(m.cache.plan, m.cache.img)
    conj.(vis).*m.cache.phases
    #return vis
end

# function _visibilities(m::ModelImage{M,I,<:NUFTCache{A}},
#                        u::AbstractArray,
#                        v::AbstractArray) where {M,I,G<:DFTAlg,A<:ObservedNUFT{G}}
#     checkuv(m.cache.alg.uv, u, v)
#     #vr = real.(m.cache.plan)
#     #vi = imag.(m.cache.plan)
#     visr = nuft(m.cache.plan.re, m.cache.img)
#     visi = nuft(m.cache.plan.im, m.cache.img)
#     vis =  StructArray{Complex{eltype(visr)}}((visr, visi))
#     #vis = m.cache.plan*m.cache.img
#     conj.(vis).*m.cache.phases
#     #return vis
# end



function _visibilities(m::ModelImage{M,I,<:NUFTCache{A}},
                      u::AbstractArray,
                      v::AbstractArray) where {M,I,A<:NUFT}
    return nocachevis(m, u, v)
end

function amplitudes(m::ModelImage{M,I,C}, u::AbstractArray, v::AbstractArray) where {M,I,C<:NUFTCache}
    vis = visibilities(m, u, v)
    return map(abs, vis)
end

function bispectra(m::ModelImage{M,I,C},
                    u1::AbstractArray,
                    v1::AbstractArray,
                    u2::AbstractArray,
                    v2::AbstractArray,
                    u3::AbstractArray,
                    v3::AbstractArray) where {M,I,C<:NUFTCache}
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
                        v3::AbstractArray) where {M,I,C<:NUFTCache}
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
                               ) where {M,I,C<:NUFTCache}
    # TODO: Fix this so we can actually use a stored cache for closures
    n = length(u1)
    u = vcat(u1,u2,u3,u4)
    v = vcat(v1,v2,v3,v4)
    vis = nocachevis(m, u, v)
    amp1 = @view vis[1:n]
    amp2 = @view vis[n+1:2n]
    amp3 = @view vis[2n+1:3n]
    amp4 = @view vis[3n+1:4n]

    lcamp = @. log(abs(amp1)*abs(amp2)/(abs(amp3)*abs(amp4)))
    return lcamp
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
Base.@kwdef struct NFFTAlg <: NUFT
    """
    Amount to pad the image
    """
    padfac::Int = 1
    """
    Controls the accuracy of the NFFT usually don't need to change this
    """
    m::Int = 10
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
