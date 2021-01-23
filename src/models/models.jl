abstract type AbstractModel{T} end
abstract type GeometricModel{T} <: AbstractModel{T} end
abstract type AbstractRadioImage{T} <: AbstractModel{T} end
abstract type AbstractPolarizedImage{T} <: AbstractModel{T} end
abstract type AbstractModifier{T} <: AbstractModel{T} end

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


function visibility(::IsNumeric, m::M, u::T, v::T, args...)

end


"""
    $(SIGNATURES)
Computes the complex visibility for the given model with args...

## Notes
THis will use FourierStyle trait to decide whether to use the numerical
Fourier transform or (if defined) the analytical Fourier transform.
"""
visibility(m::M, u::T,v::T,args...) where {M<:AbstractModel, T<:Real} =
        visibility(FourierStyle(M),m,u,v, args...)

function visibilities!(m::M,
                    vis::AbstractVector{Complex{T}},
                    u::AbstractVector{T},
                    v::AbstractVector{T},
                    args...) where {M<:AbstractModel, T<:Real}

    @inbounds for i in eachindex(u,v,vis)
        vis[i] = visibility(m, u[i], v[i], args...)
    end
    return vis
end

function visibilities(m::M,
                      u::AbstractVector{T},
                      v::AbstractVector{T},
                      args...) where {M<:AbstractModel, T<:Real}
    vis = similar(Complex{T}, u)
    return visibilities!(m, vis, u, v, args...)
end





"""
    $(SIGNATURES)
Computes the visibility amplitude of a model `m` in `Jy` and the uv positions `u`, `v`.
"""
visibility_amplitude(m::AbstractModel, u, v, args...) = abs(visibility(m, u, v, args...))


"""
    $(SIGNATURES)
Computes the complex bispectrum of a model `m` in `Jy` at the closure triangle, uv1, uv2, uv3
"""
function bispectrum(m::AbstractModel, uv1, uv2, uv3, args...)
    return visibility(m, uv1, args...)*visibility(m, uv2, args...)*visibility(m, uv3, args...)
end

"""
    $(SIGNATURES)
Computes the closure phase of a model `m` in `Jy` at closure triangle, uv1, uv2, uv3
"""
function closure_phase(m::AbstractModel, uv1, uv2, uv3, args...)
    return angle(bispectrum(m, uv1, uv2, uv3, args...))
end

"""
    $(SIGNATURES)
Computes the log closure amplitude of a model `m` in `Jy` at closure quadrangle, uv1, uv2, uv3, uv4

## Notes
We use the convention

V(uv1)*V(uv3)/V(uv2)*V(uv4),

where V is the complex visibility
"""
function logclosure_amplitude(m::AbstractModel, uv1, uv2, uv3, uv4, args...)
    return log( visibility(m, uv1...)*visibility(m, uv3...)/
               (visibility(m, uv2...)*visibility(m, uv4...))
              )
end



"""
    Returns the intensity for the model with args...
"""
function intensity(model::AbstractModel, args...) end


include("geometric_models.jl")
include("radio_image_models.jl")
include("modifiers.jl")
include("stokesmatrix.jl")
