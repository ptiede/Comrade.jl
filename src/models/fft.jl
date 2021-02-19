using FFTW
using PaddedViews
using Interpolations

struct FFTCache{T,I}
    plan::T
    xitr::I
    yitr::I
    uitr::I
    vitr::I
end

function FFTCache(sim::ROSE.StokesMatrix{T,S}) where {T,S}
    ny,nx = size(sim)
    xitr = range(-fovx/2, fovx2/2, length=nx)
    yitr = range(-fovy/2, fovy/2, length=ny)

    umax = one(T)/(2*step(xitr))
    vmax = one(T)/(2*step(yitr))

    du = umax*2/nx
    dv = vmax*2/ny

    uitr = range(-umax, umax, step=du)
    vitr = range(-vmax, vmax, step=dv)

    return FFTCache(plan_fft(sim), xitr, yitr, uitr, vitr)
end

function cft(img,x,y)
    dx = step(x)
    dy = step(y)
    umax = 1.0/(2*dx)
    vmax = 1.0/(2*dy)

    du = umax*2/(length(x))
    dv = vmax*2/(length(y))

    uu = range(-umax, umax-du, step=du)
    vv = range(-vmax, vmax-dv, step=dv)

    vis = conj.(fftshift(fft(img)))

    @inbounds @simd for I in CartesianIndices(vis)
        iy, ix = Tuple(I)
        u = -umax + du*(ix-1)
        v = -vmax + dv*(iy-1)
        vis[I] *= dx*dy*exp(2im*π*(u*first(x) + v*first(y)))
    end
    itp = interpolate(vis, BSpline(Cubic(Line(OnGrid()))))
    etp = extrapolate(itp, zero(eltype(vis)))
    sitp = Interpolations.scale(etp, uu, vv)

    return sitp
end
