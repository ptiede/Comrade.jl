abstract type AbstractModel{T} end
abstract type GeometricModel{T} <: AbstractModel{T} end
abstract type AbstractRadioImage{T} <: AbstractModel{T} end
abstract type AbstractPolarizedImage{T} <: AbstractModel{T} end

"""
    $(TYPEDEF)
Abstract type for image modifiers. These are some model wrappers
that can transform any model using simple Fourier transform properties.
To see the implemented modifier
"""
abstract type AbstractModifier{M<:AbstractModel,T} <: AbstractModel{T} end

"""
    $(SIGNATURES)
Returns a list containing all the image modifiers available in ROSE.
"""
modifierlist() = subtypes(AbstractModifier)

abstract type ImageKernel end

"""
    $(TYPEDEF)
Abstract trait for Fourier transform trait. This will decide whether
a given model should use the analytic (if it exists) or numerical
FT.
"""
abstract type FourierStyle end

struct IsAnalytic <: FourierStyle end
struct IsNumeric <: FourierStyle end


"""
    FourierStyle(::AbstractModel)
Sets the trait for the Fourier transform to be used for a given model.
The default value is IsNumeric() which says to use the numerical Fourier
transform.
"""
FourierStyle(::Type{<:AbstractModel}) = IsNumeric()

FourierStyle(::Type{<:GeometricModel}) = IsAnalytic()

FourierStyle(::Type{<:AbstractModifier}) = IsAnalytic()

abstract type ObservationCache end
abstract type ModelCache <: ObservationCache end

struct NFFTCache{T} <: ModelCache
    cache::T
end
struct FFTCache{T} <: ModelCache
    cache::T
end
struct DFTCache{T} <: ModelCache
    cache::T
end

abstract type CacheStyle end

struct WithCache <: CacheStyle end
struct NoCache <: CacheStyle end

CacheStyle(::Type{<:AbstractModel}) = NoCache();
CacheStyle(::AbstractModifier{T,M}) where {T,M} = CacheStyle(M)

"""
    $(SIGNATURES)
Creates a model cache to speed up computation. For certain models, e.g.
numerical visibility models with fix u,v, positions, e.g. models that don't
change the uv positions for certain parameters, this can cause a speed up.
"""
function createcache(m::M, obs::Observation, cachetype::ObservationCache) where {M<: AbstractModel}
    return createcache(CacheStyle(M), m, obs, cachetype)
end




"""
    $(SIGNATURES)
Computes the complex visibility for the given model with args...

## Notes
THis will use FourierStyle trait to decide whether to use the numerical
Fourier transform or (if defined) the analytical Fourier transform.
"""
visibility(m::M, u,v,args...) where {M<:AbstractModel} =
        visibility(FourierStyle(M),m,u,v, args...)


function visibility(::IsNumeric, m::M, u, v, args...) where {T,M<:AbstractModel{T}}
    throw("Not implemented yet")
end





function visibilities!(
                    vis::AbstractVector{Complex{S}},
                    m::M,
                    u::AbstractVector,
                    v::AbstractVector,
                    args...) where {S,M<:AbstractModel}

    @inbounds for i in eachindex(u,v,vis)
        vis[i] = visibility(m, u[i], v[i], args...)
    end
    return vis
end

function visibilities(m::M,
                      u::AbstractVector,
                      v::AbstractVector,
                      args...) where {S,M<:AbstractModel{S}}
    vis = StructArray{Complex{S}}(re=similar(u,S), im=similar(v,S))
    return visibilities!(vis, m, u, v, args...)
end





"""
    $(SIGNATURES)
Computes the visibility amplitude of a model `m` in `Jy` and the uv positions `u`, `v`.
"""
@inline visibility_amplitude(m::AbstractModel, u, v, args...) = abs(visibility(m, u, v, args...))


"""
    $(SIGNATURES)
Computes the complex bispectrum of a model `m` in `Jy` at the closure triangle, uv1, uv2, uv3
"""
@inline function bispectrum(m::AbstractModel, u1,v1, u2, v2, u3, v3, args...)
    return visibility(m, u1, v1, args...)*
           visibility(m, u2, v2, args...)*
           visibility(m, u3, v3, args...)
end

"""
    $(SIGNATURES)
Computes the closure phase of a model `m` in `Jy` at closure triangle, uv1, uv2, uv3
"""
@inline function closure_phase(m::AbstractModel, u1, v1, u2, v2, u3, v3, args...)
    return angle(bispectrum(m, u1,v1, u2,v2, u3,v3, args...))
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


function intensitymap!(im::StokesImage{T,S}, m::AbstractModel) where {T,S}
    ny,nx = size(im)
    psizex = im.fovx/max(nx-1,1)
    psizey = im.fovy/max(ny-1,1)
    flux = zero(T)
    fov = max(im.fovx, im.fovy)
    npix = max(nx, ny)
    @inbounds @simd for I in CartesianIndices(im)
        iy,ix = Tuple(I)
        x = -im.fovx/2 + psizex*(ix-1)
        y = -im.fovy/2 + psizey*(iy-1)
        tmp = intensity(m, x, y, fov, npix)
        im[I] = tmp
    end
    return im
end

function intensitymap(m::AbstractModel{T}, nx, ny, fovx, fovy) where {T}
    im = zeros(T, ny, nx)
    sim = StokesImage(im, fovx, fovy)
    intensitymap!(sim, m)

end





"""
    intensity(model::AbstractModel, x, y, args...)
Returns the intensity for the model at x and y with args...
"""
function intensity end

include("geometric_models.jl")
include("radio_image_models.jl")
include("modifiers.jl")
