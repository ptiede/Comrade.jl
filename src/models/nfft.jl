function create_plan(img, u, v, padfac, m)
    dx,dy = pixelsizes(img)
    uv = Matrix{eltype(u)}(undef, 2, length(u))
    uv[1,:] = u*dx
    uv[2,:] = v*dy
    dx,dy = pixelsizes(img)
    ny,nx = size(img)
    nnx = nextpow(2, padfac*nx)
    nny = nextpow(2, padfac*ny)
    nsx = nnx÷2-nx÷2
    nsy = nny÷2-ny÷2
    pimg = PaddedView(zero(eltype(img)), img.im,
                      (1:nnx, 1:nny),
                      (nsx+1:nsx+nx, nsy+1:nsy+ny)
                     )
    plan = plan_nfft(uv, size(pimg), m)
    phases = @. exp(1im*π*(u*dx + v*dy))*dx*dy#*exp(1im*π*(u*x + v*y))*dx*dy
    return plan, phases, pimg
end

Base.@kwdef struct NFFT <: FourierTransform
    padfac::Int = 2
    m::Int = 10
end

struct NFFTCache{C} <: AbstractCache
    plan::C
    padfac::Int
    m::Int
end
