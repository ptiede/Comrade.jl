using Comrade
using ForwardDiff
using AbstractFFTs
using FFTW
using PaddedViews
# @edit fft(x_plus_dx, 1:1) points me here:
#AbstractFFTs.complexfloat(x::AbstractArray{<:ForwardDiff.Dual}) = float.(ForwardDiff.value.(x) .+ 0im)

# @edit fft(x_plus_dx .+ 0im, 1:1) # now this makes a plan, we need:
#AbstractFFTs.plan_fft(x::AbstractArray{<:ForwardDiff.Dual}, region=1:ndims(x)) = plan_fft(ForwardDiff.value.(x) .+ 0im, region)
#AbstractFFTs.plan_fft(x::AbstractArray{<:Complex{<:ForwardDiff.Dual}}, region=1:ndims(x)) = plan_fft(ForwardDiff.value.(x), region)

# Where I want value() to work on complex duals too:
#ForwardDiff.value(x::Complex{<:ForwardDiff.Dual}) = Complex(x.re.value, x.im.value)
#ForwardDiff.partials(x::Complex{<:ForwardDiff.Dual}, n::Int) = Complex(ForwardDiff.partials(x.re, n), ForwardDiff.partials(x.im, n))
#ForwardDiff.npartials(x::Complex{<:ForwardDiff.Dual}) = ForwardDiff.npartials(x.re)

# Now fft(x_plus_dx) fails at *(p::FFTW.cFFTWPlan{ComplexF64, -1, false, 1, UnitRange{Int64}}, x::Vector{Complex{ForwardDiff.Dual{Nothing, Float64, 2}}}), great!

function Base.:*(p::AbstractFFTs.Plan, x::PaddedView{<:ForwardDiff.Dual{T,V,P},N, I,<:IntensityMap}) where {T,V,P,N,I}
    M = typeof(ForwardDiff.value(first(x)))
    cache = Matrix{M}(undef, size(x)...)
    cache .= ForwardDiff.value.(x)
    xtil = p * cache
    ndxs = ForwardDiff.npartials(first(x))
    dxtils = ntuple(ndxs) do n
        cache[:,:] .= ForwardDiff.partials.(x, n)
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



function f(r,α)
    m = ExtendedRing(r, α)
    img = intensitymap(m, 100.0, 100.0, 128,128)
    dx,dy = pixelsizes(img)
    ny,nx = size(img)
    padfac = 2
    nnx = nextpow(2, padfac*nx)
    nny = nextpow(2, padfac*ny)
    pimg = PaddedView(zero(eltype(img)), img, (nnx, nny))
    uu = fftshift(fftfreq(nnx, 1/dx))
    vv = fftshift(fftfreq(nny, 1/dy))
    #println(img.im[1])
    #plan = plan_fft(pimg)
    #println(img.im[1])
    #println(typeof(pimg))
    p = plan_fft(pimg)
    vis = fftshift(p*pimg)

    return sum(abs2, vis)
end

#x = rand(10)
#
#xp = plan_fft(x)
#
#xhat = xp*x
#
#x_plus_dx = [ForwardDiff.Dual(x[i], (i,i^2)) for i in 1:10]
#dx1 = ForwardDiff.partials.(x_plus_dx, 1)
#dx1hat = xp * dx1
#dx2hat = xp * @view reinterpret(Float64, x_plus_dx)[3:3:end]
#xtil_plus = [Complex(   # re-assemble
#             ForwardDiff.Dual(real(xhat[i]), (real(dx1hat[i]), real(dx2hat[i]))),
#             ForwardDiff.Dual(imag(xhat[i]), (imag(dx1hat[i]), imag(dx2hat[i])))
#            ) for i in 1:10][1]
#
#fft(x_plus_dx)[1]
