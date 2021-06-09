struct ModelCache{P,M} <: ObservationCache
    cache::P
    model::M
end


struct FFTPlan{T,S,I,V,C}
    plan::T
    sim::S
    uitr::I
    vitr::I
    u::V
    v::V
    phases::C
end

function FFTPlan(sim::ROSE.StokesImage{T,S}, u, v) where {T,S}
    ny,nx = size(sim)
    x,y = imagepixels(sim)
    dx,dy = pixelsizes(sim)

    uu = fftshift(fftfreq(length(x), 1/dx))
    vv = fftshift(fftfreq(length(y), 1/dy))

    return FFTCache(plan_fft(sim), x, y, uu, vv, u, v)
end


struct NFFTPlan{T,P,V,C} <: FTCache
    plan::T
    sim::P
    u::V
    v::V
    phases::C
end

function NFFTPlan(image::StokesImage, u, v)

    @assert length(u) == length(v) "u and v must have equal length"

    # Get pixels sizes since needed to switch to grid in NFFT
    dx, dy = ROSE.pixelsizes(image)
    ny, nx = size(image)


    uv = zeros(2, length(u))
    uv[1,:] = u*dx
    uv[2,:] = v*dy

    p = plan_nfft(uv, (ny, nx))
    phases = exp.(1im*π*(u*dx + v*dy))*dx*dy
    NFFTCache(p, deepcopy(sim), phases)
end


struct NFFTV end
struct FFTWV end

"""
    simulate_vis
Computes the visibilities of a stokes image along the vector u and v and using an
algorithm alg
"""
function simulate_vis end

function simulate_vis(m::AbstractModel, cache::FTCache)
    image = intensitymap!(cache.sim, m)
    return simulate_vis(image, cache)
end


function simulate_vis(sim::StokesImage, u::AbstractVector, v::AbstractVector, ::NFFTV)
    dx, dy = pixelsizes(sim)
    ny, nx = size(sim)


    uv = zeros(2, length(u))
    uv[1,:] = u*dx
    uv[2,:] = v*dy

    phases = exp.(1im*π*(u*dx + v*dy))*dx*dy
    vis = nfft(uv, sim).*phases
    return vis
end

function simulate_vis(sim::StokesImage, cache::NFFTCache)
    vis =  nfft(cache.plan, sim)
    return conj.(vis).*cache.phases
end

function simulate_vis(sim::StokesImage, cache::FFTCache)
    x,y = imagepixels(sim)
    dx,dy = pixelsizes(sim)

    uu = cache.uitr
    vv = cache.vitr


    vis = fftshift(fft(sim))

    @avx for I in CartesianIndices(vis)
        iy, ix = Tuple(I)
        vis[I] =  conj(vis[I])*dx*dy*exp(2im*π*(uu[ix]*first(x) + vv[iy]*first(y)))
    end
    itp = interpolate(vis, BSpline(Cubic(Line(OnGrid()))))
    etp = extrapolate(itp, zero(eltype(vis)))
    sitp = scale(etp, cache.uitr, cache.vitr)

    return sitp.(cache.u,cache.v)
end

function simulate_vis(sim::StokesImage, u::AbstractVector, v::AbstractVector, ::FFTWV)
    dx,dy = pixelsizes(sim)
    x,y = imagepixels(sim)
    fu = fftshift(fftfreq(length(x), 1/dx))
    fv = fftshift(fftfreq(length(y), 1/dy))

    vis = fftshift(fft(sim))

    @inbounds for I in CartesianIndices(vis)
        iy, ix = Tuple(I)
        vis[I] =  conj(vis[I])*dx*dy*exp(2im*π*(fu[ix]*first(x) + fv[iy]*first(y)))
    end
    itp = interpolate(vis, BSpline(Cubic(Line(OnGrid()))))
    etp = extrapolate(itp, zero(eltype(vis)))
    sitp = scale(etp, fu, fv)

    return sitp.(u,v)
end
