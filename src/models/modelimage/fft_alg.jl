export FFTAlg






#function ChainRulesCore.rrule(::typeof(create_interpolator), u, v, vis)
#    y = create_interpolator(u, v, vis)
#    function create_interp_pullback(Δy)
#        Δf = NoTangent()
#        Δu = ZeroTanget()
#        Δv = ZeroTangent()
#        Δvis = Δy*create_interpolator(u, v, vis)
#        return Δf, Δu, Δv, Δvis
#    end
#    return y, create_interp_pullback
#end
#
#@adjoint create_interpolator(u, v, vis) = create_interpolator(u, v, vis), Δy->(u, v, Δy)

"""
    $(TYPEDEF)
The cache used when the `FFT` algorithm is used to compute
visibilties. This is an internal type and is not part of the public API
"""
struct FFTCache{A<:FFTAlg,P,I,S} <: AbstractCache
    alg::A # FFT algorithm
    plan::P # FFT plan or matrix
    img::I # image buffer used for FT
    sitp::S # interpolator function used to find vis at abritrary uv position.
end

# These are all overloads to allow us to forward propogate ForwardDiff dual numbers through
# an abstract FFT.
AbstractFFTs.complexfloat(x::AbstractArray{<:ForwardDiff.Dual}) = float.(ForwardDiff.value.(x) .+ 0.0im)

# Make a plan with Dual numbers
AbstractFFTs.plan_fft(x::AbstractArray{<:ForwardDiff.Dual}, region=1:ndims(x)) = plan_fft(zeros(ComplexF64, size(x)), region)
AbstractFFTs.plan_fft(x::AbstractArray{<:Complex{<:ForwardDiff.Dual}}, region=1:ndims(x)) = plan_fft(zeros(ComplexF64, size(x)), region)

# Allow things to work with Complex Dual numbers.
ForwardDiff.value(x::Complex{<:ForwardDiff.Dual}) = Complex(x.re.value, x.im.value)
ForwardDiff.partials(x::Complex{<:ForwardDiff.Dual}, n::Int) = Complex(ForwardDiff.partials(x.re, n), ForwardDiff.partials(x.img, n))
ForwardDiff.npartials(x::Complex{<:ForwardDiff.Dual}) = ForwardDiff.npartials(x.re)

# internal function that creates the interpolator objector to evaluate the FT.
function create_interpolator(U, V, vis::AbstractArray{<:Complex}, pulse)
    # Construct the interpolator
    #itp = interpolate(vis, BSpline(Cubic(Line(OnGrid()))))
    #etp = extrapolate(itp, zero(eltype(vis)))
    #scale(etp, u, v)
    p1 = BicubicInterpolator(U, V, real(vis), NoBoundaries())
    p2 = BicubicInterpolator(U, V, imag(vis), NoBoundaries())
    function (u,v)
        pl = visibility_point(pulse, u, v, 0.0, 0.0)
        return pl*(p1(u,v) + 1im*p2(u,v))
    end
end

function create_interpolator(U, V, vis::StructArray{<:StokesParams}, pulse)
    # Construct the interpolator
    pI_real = BicubicInterpolator(U, V, real(vis.I), NoBoundaries())
    pI_imag = BicubicInterpolator(U, V, real(vis.I), NoBoundaries())

    pQ_real = BicubicInterpolator(U, V, real(vis.Q), NoBoundaries())
    pQ_imag = BicubicInterpolator(U, V, real(vis.Q), NoBoundaries())

    pU_real = BicubicInterpolator(U, V, real(vis.U), NoBoundaries())
    pU_imag = BicubicInterpolator(U, V, real(vis.U), NoBoundaries())

    pV_real = BicubicInterpolator(U, V, real(vis.V), NoBoundaries())
    pV_imag = BicubicInterpolator(U, V, real(vis.V), NoBoundaries())


    function (u,v)
        pl = visibility_point(pulse, u, v, 0.0, 0.0)
        return StokesParams(
            pI_real(u,v)*pl + 1im*pI_imag(u,v)*pl,
            pQ_real(u,v)*pl + 1im*pQ_imag(u,v)*pl,
            pU_real(u,v)*pl + 1im*pU_imag(u,v)*pl,
            pV_real(u,v)*pl + 1im*pV_imag(u,v)*pl,
        )
    end
end

#=
This is so ForwardDiff works with FFTW. I am very harsh on the `x` type because of type piracy.
=#
function Base.:*(p::AbstractFFTs.Plan, x::PaddedView{<:ForwardDiff.Dual{T,V,P},N, I,<:IntensityMapTypes}) where {T,V,P,N,I}
    M = typeof(ForwardDiff.value(first(x)))
    cache = Matrix{M}(undef, size(x)...)
    cache .= ForwardDiff.value.(x)
    xtil = p * cache
    ndxs = ForwardDiff.npartials(first(x))
    dxtils = ntuple(ndxs) do n
        cache .= ForwardDiff.partials.(x, n)
        p * cache
    end
    out = similar(cache, Complex{ForwardDiff.Dual{T,V,P}})
    for i in eachindex(out)
        dual = getindex.(dxtils, i)
        prim = xtil[i]
        red = ForwardDiff.Dual{T,V,P}(real(prim), ForwardDiff.Partials(real.(dual)))
        imd = ForwardDiff.Dual{T,V,P}(imag(prim), ForwardDiff.Partials(imag.(dual)))
        out[i] = Complex(red, imd)
    end
    return out
end

# function padimage(alg::FFTAlg, img)
#     padfac = alg.padfac
#     ny,nx = size(img)
#     nnx = nextpow(2, padfac*nx)
#     nny = nextpow(2, padfac*ny)
#     nsx = nnx÷2-nx÷2
#     nsy = nny÷2-ny÷2
#     return PaddedView(zero(eltype(img)), img,
#                       (1:nnx, 1:nny),
#                       (nsx+1:nsx+nx, nsy+1:nsy+ny)
#                      )
# end

function padimage(img::IntensityMap, alg::FFTAlg)
    padfac = alg.padfac
    ny,nx = size(img)
    nnx = nextprod((2,3,5,7), padfac*nx)
    nny = nextprod((2,3,5,7), padfac*ny)
    PaddedView(zero(eltype(img)), img, (nny, nnx))
end

function padimage(img::StokesIntensityMap, alg::FFTAlg)
    pI = padimage(stokes(img, :I), alg)
    pQ = padimage(stokes(img, :Q), alg)
    pU = padimage(stokes(img, :U), alg)
    pV = padimage(stokes(img, :V), alg)
    return StructArray{eltype(img)}((I=pI, Q=pQ, U=pU, V=pV))
end

# phasecenter the FFT.
function phasecenter(vis, X, Y, U, V)
    x0 = first(X)
    y0 = first(Y)
    return conj.(vis).*cispi.(2 .* (U.*x0 .+ V'.*y0))
end

function applyfft(plan, img::AbstractArray{<:Number})
    return fftshift(plan*img)
end

function applyfft(plan, img::AbstractArray{<:StokesParams})
    visI = applyfft(plan, stokes(img, :I))
    visQ = applyfft(plan, stokes(img, :Q))
    visU = applyfft(plan, stokes(img, :U))
    visV = applyfft(plan, stokes(img, :V))
    return StructArray{StokesParams{eltype(visI)}}((I=visI, Q=visQ, U=visU, V=visV))
end



function create_cache(alg::FFTAlg, img::IntensityMapTypes, pulse::Pulse=DeltaPulse())
    pimg = padimage(img, alg)
    # Do the plan and then fft
    plan = plan_fft(pimg)
    vis = applyfft(plan, pimg)

    #Construct the uv grid
    (;X, Y) = imagepixels(img)
    (;U, V) = uviterator(size(pimg, 1), step(X), size(pimg, 2), step(Y))


    vispc = phasecenter(vis, X, Y, U, V)
    sitp = create_interpolator(U, V, vispc, stretched(pulse, step(X), step(Y)))
    return FFTCache(alg, plan, img, sitp)
end

function update_cache(cache::FFTCache, img::IntensityMapTypes, pulse)
    plan = cache.plan
    pimg = padimage(img, cache.alg)
    vis = applyfft(plan, pimg)

    (;X,Y) = imagepixels(img)

    # Flip U because of radio conventions
    (;U, V) = uviterator(size(pimg, 1), step(X), size(pimg, 2), step(Y))

    vispc = phasecenter(vis, X, Y, U, V)
    sitp = create_interpolator(uu, vv, vispc, img)
    return FFTCache(cache.alg, plan, img, sitp)
end


"""
    uviterator(nx, dx, ny dy)

Construct the u,v iterators for the Fourier transform of the image
with pixel sizes `dx, dy` and number of pixels `nx, ny`

If you are extending Fourier transform stuff please use these functions
to ensure that the centroid is being properly computed.
"""
function uviterator(nx, dx, ny, dy)
    U = fftshift(fftfreq(nx, inv(dx)))
    V = fftshift(fftfreq(ny, inv(dy)))
    return (;U, V)
end

#function ChainRulesCore.rrule(::typeof(phasecenter!), vis, uu, vv, x0, y0, dx, dy)
#    vis = phasecenter!(vis, uu, vv, x0, y0, dx, dy)
#    function phasecenter!_pullback(Δy)
#        Δf = NoTangent()
#        Δvis = @thunk(dx*dy*exp.(2im*π*(uu*x0 + vv'*y0)))
#        Δu = @thunk(dx*dy*2im*π*x0*exp.(2im*π*(uu*x0 + vv'*y0)))
#        Δv = @thunk(dx*dy*2im*π*y0*exp.(2im*π*(uu*x0 + vv'*y0)))
#    end
#end

"""
    fouriermap(m, x)

Create a Fourier or visibility map of a model `m`
where the image is specified in the image domain by the
pixel locations `x` and `y`
"""
function fouriermap(m, dims::AbstractDims)
    X = dims.X
    Y = dims.Y
    uu,vv = uviterator(length(X), step(X), length(Y), step(Y))
    uvgrid = ComradeBase.grid(U=uu, V=vv)
    vis = visibility.(Ref(m), uvgrid)

    return vis
end




# function fouriermap(m::ModelImage, fovx, fovy, x0, y0, nx, ny)
#     cache = create_cache(FFTAlg(), m.image)
#     x,y = imagepixels(fovx, fovy, x0, y0, nx, ny)
#     dx = step(x); dy = step(y)
#     uu,vv = uviterator(dx, dy, nx, ny)

#     T = Complex{eltype(m.image)}
#     vis = Matrix{T}(undef, ny, nx)

#     @inbounds for I in CartesianIndices(vis)
#         iy, ix = Tuple(I)
#         vp = cache.sitp(uu[ix], vv[iy])
#         vis[I] = vp
#     end
#     return vis

# end

function phasedecenter!(vis, X, Y)
    nx = length(X)
    ny = length(Y)
    uu,vv = uviterator(nx, step(X), ny, step(Y))
    x0 = first(X)
    y0 = first(Y)
    for I in CartesianIndices(vis)
        iy, ix = Tuple(I)
        vis[I] = conj(vis[I])*cispi(2*(uu[ix]*x0 + vv[iy]*y0))*nx*ny
    end
    return vis
end

# function _visibilities(mimg::ModelImage{M, I, <:FFTCache}, u, v, args...) where {M,I}
#     return visibility.(Ref(mimg), u, v, args...)
# end

@inline function visibility_point(mimg::ModelImage{M,I,<:FFTCache}, u, v, time, freq) where {M,I}
    return mimg.cache.sitp(u, v)
end
