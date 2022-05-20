export NFFTAlg

using NFFT

function NFFTAlg(obs::EHTObservation; kwargs...)
    u, v = getuv(arrayconfig(obs))
    return NFFTAlg(u, v; kwargs...)
end

function NFFTAlg(u::AbstractArray, v::AbstractArray; padfac=1, m=10)
    uv = Matrix{eltype(u)}(undef, 2, length(u))
    uv[1,:] .= u
    uv[2,:] .= v
    return ObservedNUFT(NFFTAlg(;padfac, m), uv)
end

function NFFTAlg(ac::ArrayConfiguration; padfac=1, m=10)
    u, v = getuv(ac)
    return NFFTAlg(u, v; padfac, m)
end


function padimage(alg::NFFTAlg, img)
    padfac = alg.padfac
    # if no padding exit now
    (padfac == 1) && return img

    ny,nx = size(img)
    nnx = nextpow(2, padfac*nx)
    nny = nextpow(2, padfac*ny)
    nsx = nnx÷2-nx÷2
    nsy = nny÷2-ny÷2
    return PaddedView(zero(eltype(img)), img,
                      (1:nnx, 1:nny),
                      (nsx+1:nsx+nx, nsy+1:nsy+ny)
                     )
end

function plan_nuft(alg::ObservedNUFT{<:NFFTAlg}, img, dx, dy)
    uv2 = similar(alg.uv)
    uv2[1,:] .= alg.uv[1,:]*dx
    uv2[2,:] .= alg.uv[2,:]*dy
    plan = plan_nfft(uv2, size(img); precompute=NFFT.POLYNOMIAL)
    return plan
end

function make_phases(alg::ObservedNUFT{<:NFFTAlg}, img)
    dx, dy = pixelsizes(img)
    u = @view alg.uv[1,:]
    v = @view alg.uv[2,:]
    return cispi.((u.*dx .+ v.*dy)).*visibility_point.(Ref(img.pulse), u.*dx, v.*dy)
end

@inline function create_cache(alg::ObservedNUFT{<:NFFTAlg}, plan, phases, img)
    timg = IntensityMap(transpose(img.im), img.fovx, img.fovy, img.pulse)
    return NUFTCache(alg, plan, phases, img.pulse, timg)
end

function _frule_vis(m::ModelImage{M,<:IntensityMap{<:ForwardDiff.Dual{T,V,P}},<:NUFTCache{O}}) where {M,T,V,P,A<:NFFTAlg,O<:ObservedNUFT{A}}
    S = typeof(ForwardDiff.value(first(m.cache.img)))
    p = m.cache.plan
    # Compute the fft
    buffer = similar(m.cache.img, S)
    buffer .= ForwardDiff.value.(m.cache.img)
    xtil = p*buffer
    out = similar(buffer, Complex{ForwardDiff.Dual{T,V,P}})
    # Now take the deriv of nuft
    ndxs = ForwardDiff.npartials(first(m.cache.img))
    dxtils = ntuple(ndxs) do n
        buffer .= ForwardDiff.partials.(m.cache.img, n)
        p * buffer
    end
    out = similar(xtil, Complex{ForwardDiff.Dual{T,V,P}})
    for i in eachindex(out)
        dual = getindex.(dxtils, i)
        prim = xtil[i]
        red = ForwardDiff.Dual{T,V,P}(real(prim), ForwardDiff.Partials(real.(dual)))
        imd = ForwardDiff.Dual{T,V,P}(imag(prim), ForwardDiff.Partials(imag.(dual)))
        out[i] = Complex(red, imd)
    end
    return out
end

function _visibilities(m::ModelImage{M,<:IntensityMap{<:ForwardDiff.Dual{T,V,P}},<:NUFTCache{O}},
    u::AbstractArray,
    v::AbstractArray) where {M,T,V,P,A<:NFFTAlg,O<:ObservedNUFT{A}}
    checkuv(m.cache.alg.uv, u, v)
    # Now reconstruct everything

    vis = _frule_vis(m)
    conj.(vis).*m.cache.phases
#return vis
end
