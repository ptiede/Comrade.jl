export NFFTAlg

using NFFT

"""
    NFFTAlg(obs::EHTObservation; padfac=1, m=10)

Create an algorithm object using the non-unform Fourier transform object from the observation
`obs`. This will extract the uv positions from the observation to allow for a more efficient
FT cache.

The optional arguments are: `padfac` specifies how much to pad the image by, and `m`
is an internal variable for `NFFT.jl`.
"""
function NFFTAlg(obs::EHTObservation; padfac=1, m=10)
    u, v = getuv(arrayconfig(obs))
    return NFFTAlg(u, v; padfac, m)
end

"""
    NFFTAlg(u::AbstractArray, v::AbstractArray; padfac=1, m=10)

Create an algorithm object using the non-unform Fourier transform object from uv positions
`u`, `v`. This will extract the uv positions from the observation to allow for a more efficient
FT cache.

The optional arguments are: `padfac` specifies how much to pad the image by, and `m`
is an internal variable for `NFFT.jl`.
"""
function NFFTAlg(u::AbstractArray, v::AbstractArray; padfac=1, m=10)
    uv = Matrix{eltype(u)}(undef, 2, length(u))
    uv[1,:] .= u
    uv[2,:] .= v
    return ObservedNUFT(NFFTAlg(;padfac, m), uv)
end

"""
    NFFTAlg(ac::ArrayConfiguration; padfac=1, m=10)

Create an algorithm object using the non-unform Fourier transform object from the array
configuration `ac`. This will extract the uv positions from the observation to allow
for a more efficient FT cache.

The optional arguments are: `padfac` specifies how much to pad the image by, and `m`
is an internal variable for `NFFT.jl`.
"""
function NFFTAlg(ac::ArrayConfiguration; padfac=1, m=10)
    u, v = getuv(ac)
    return NFFTAlg(u, v; padfac, m)
end


# pad from the center of the position.
function padimage(alg::NFFTAlg, img)
    padfac = alg.padfac
    # if no padding exit now
    (padfac == 1) && return img

    ny,nx = size(img)
    nnx = nextprod((2,3,5,7), padfac*nx)
    nny = nextprod((2,3,5,7), padfac*ny)
    nsx = nnx÷2-nx÷2
    nsy = nny÷2-ny÷2
    pimg = PaddedView(zero(eltype(img)), img.im,
                      (1:nnx, 1:nny),
                      (nsx+1:nsx+nx, nsy+1:nsy+ny)
                     )
    dx, dy = pixelsizes(img)
    return IntensityMap(collect(pimg), dx*size(pimg,2), dy*size(pimg, 1), img.pulse)
end

function plan_nuft(alg::ObservedNUFT{<:NFFTAlg}, img, dx, dy)
    uv2 = similar(alg.uv)
    uv2[1,:] .= alg.uv[1,:]*dx
    uv2[2,:] .= alg.uv[2,:]*dy
    plan = plan_nfft(uv2, size(img'); precompute=alg.alg.precompute)
    return plan
end

function make_phases(alg::ObservedNUFT{<:NFFTAlg}, img)
    dx, dy = pixelsizes(img)
    u = @view alg.uv[1,:]
    v = @view alg.uv[2,:]
    return cispi.((u.*dx .+ v.*dy)).*visibility_point.(Ref(img.pulse), u.*dx, v.*dy)
end

@inline function create_cache(alg::ObservedNUFT{<:NFFTAlg}, plan, phases, img)
    #timg = #IntensityMap(transpose(img.im), img.fovx, img.fovy, img.pulse)
    return NUFTCache(alg, plan, phases, img.pulse, img.im')
end

# Allow NFFT to work with ForwardDiff.
function _frule_vis(m::ModelImage{M,<:IntensityMap{<:ForwardDiff.Dual{T,V,P}},<:NUFTCache{O}}) where {M,T,V,P,A<:NFFTAlg,O<:ObservedNUFT{A}}
    S = typeof(ForwardDiff.value(first(m.cache.img)))
    p = m.cache.plan
    # Compute the fft
    buffer = similar(m.cache.img, S)
    buffer .= ForwardDiff.value.(m.cache.img)
    xtil = p*(buffer .+ 0.0im)
    out = similar(buffer, Complex{ForwardDiff.Dual{T,V,P}})
    # Now take the deriv of nuft
    ndxs = ForwardDiff.npartials(first(m.cache.img))
    dxtils = ntuple(ndxs) do n
        buffer .= ForwardDiff.partials.(m.cache.img, n)
        p * (buffer .+ 0.0im)
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
    return conj.(vis).*m.cache.phases
end
