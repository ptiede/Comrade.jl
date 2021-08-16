abstract type AbstractModel{T} end
abstract type ObservationCache end

include("traits.jl")




function intensity(m::M, u, v, args...) where {M<:AbstractModel}
    return intensity(ImStyle(M), m, u, v)
end


function intensitymap!(im::StokesImage{T,S}, m::AbstractModel) where {T,S}
    ny,nx = size(im)
    psizex = im.fovx/max(nx-1,1)
    psizey = im.fovy/max(ny-1,1)
    fov = max(im.fovx, im.fovy)
    npix = max(nx, ny)
    @inbounds @simd for I in CartesianIndices(im)
        iy,ix = Tuple(I)
        x = -im.fovx/2 + psizex*(ix-1)
        y = -im.fovy/2 + psizey*(iy-1)
        tmp = intensity(m, x, y, im.fovx, im.fovy, nx, ny)
        im[I] = tmp
    end
    return im
end


@memoize function imcache(m::AbstractModel, fovx, fovy, nx, ny)
    sim = ROSE.StokesImage(zeros(ny,nx), fovx, fovy)
    sim = intensitymap!(sim, m)
    x_itr,y_itr = imagepixels(sim)
    itp = interpolate(sim, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, x_itr, y_itr)
    etp = extrapolate(sitp, 0)
    return etp
end

function intensity(::ImNumeric, m::AbstractModel, x, y, fovx=120.0,fovy=120.0,nx=200,ny=200, args...)
    itp = imcache(m, fovx, fovy, nx, ny)
    return itp(x,y)
end


function intensitymap(m::AbstractModel{T}, nx, ny, fovx, fovy) where {T}
    im = zeros(T, ny, nx)
    sim = StokesImage(im, fovx, fovy)
    intensitymap!(sim, m)
end




function visibilities(::VisNumeric, m::AbstractModel, cache::FTCache)
    return simulate_vis(m, cache)
end

function visibilities(::VisNumeric, m::AbstractModel,
                      u::AbstractVector, v::AbstractVector,
                      alg=NFFTV();
                      fov=200.0, npix = 512)

    sim = intensitymap(m, npix, npix, fov, fov)

end

function visibilities(m::AbstractModel, cache, args...)
    return visibilities(VisStyle(M), m, cache, args...)
end


function visibilities(m::AbstractModel,
                      u::AbstractVector,
                      v::AbstractVector,
                      args...)
        return visibilities(VisStyle(M), m, u, v, args...)
end

function visibilities(::VisAnalytic, m::AbstractModel, u, v, args...)
    return mappedarray(eachindex(u,v)) do i
        vis[i] = visibility(m, u[i], v[i], args...)
    end
end




"""
    $(SIGNATURES)
Computes the visibility amplitude of a model `m` in `Jy` and the uv positions `u`, `v`.
"""
@inline visibility_amplitudes(m::AbstractModel, u, v, args...) = abs.(visibilities(m, u, v, args...))


"""
    $(SIGNATURES)
Computes the complex bispectrum of a model `m` in `Jy` at the closure triangle, uv1, uv2, uv3
"""
@inline function bispectra(m::AbstractModel, u1,v1, u2, v2, u3, v3, args...)
    return visibilities(m, u1, v1, args...).*
           visibilities(m, u2, v2, args...).*
           visibilities(m, u3, v3, args...)
end

"""
    $(SIGNATURES)
Computes the closure phase of a model `m` in `Jy` at closure triangle, uv1, uv2, uv3
"""
@inline function closure_phases(m::AbstractModel, u1, v1, u2, v2, u3, v3, args...)
    return mappedarray(u1,v1,u2,v2,u3,v3, args...) do pu1, pv1, pu2, pv2, pu3, pv3
        angle(bispectra(m, u1[i],v1[i], u2[i] ,v2[i] , u3,v3, args...))
end

"""
    $(SIGNATURES)
Computes the log closure amplitude of a model `m` in `Jy` at closure quadrangle, uv1, uv2, uv3, uv4

## Notes
We use the convention

V(uv1)*V(uv3)/V(uv2)*V(uv4),

where V is the complex visibility
"""
@inline function logclosure_amplitude(m::AbstractModel,
                                      u1, v1,
                                      u2, v2,
                                      u3, v3,
                                      u4, v4,
                                      args...)
    return log( visibility_amplitude(m, u1,v1...)*visibility_amplitude(m, u2,v3...)/
               (visibility_amplitude(m, u3,v2...)*visibility_amplitude(m, u4,v4...))
              )
end









"""
    intensity(model::AbstractModel, x, y, args...)
Returns the intensity for the model at x and y with args...
"""
function intensity end

include("geometric_models.jl")
include("radio_image_models.jl")
include("composite.jl")


include("modifiers.jl")
