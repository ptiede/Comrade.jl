export NFFTAlg

using NFFT
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


function padimage(img, alg::NFFTAlg)
    padfac = alg.padfac
    # if no padding exit now
    (padfac == 1) && return img

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

function plan_nuft(alg::ObservedNUFT{<:NFFTAlg}, img, dx, dy)
    uv2 = copy(alg.uv)
    uv2[1,:] .= uv2[1,:]*dx
    uv2[2,:] .= uv2[2,:]*dy
    plan = plan_nfft(uv2, size(img); precompute=NFFT.POLYNOMIAL)
    return plan
end

function make_phases(alg::ObservedNUFT{<:NFFTAlg}, img)
    dx, dy = pixelsizes(img)
    u = @view alg.uv[1,:]
    v = @view alg.uv[2,:]
    return cispi.(-(u.*dx .+ v.*dy)).*visibility_point.(Ref(img.pulse), u.*dx, v.*dy).*dx.*dy
end

@inline function create_cache(alg::ObservedNUFT{<:NFFTAlg}, plan, phases, img)
    return NUFTCache(alg, plan, phases, img)
end
