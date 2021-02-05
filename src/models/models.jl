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
abstract type AbstractModifier{T} <: AbstractModel{T} end

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


function visibility(::IsNumeric, m::M, u::T, v::T, args...) where {T,M<:AbstractModel{T}}
    throw("Numeric visibility is not implemented yet")
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


function visibilities!(
                    vis::AbstractVector{Complex{S}},
                    m::M,
                    u::AbstractVector{T},
                    v::AbstractVector{T},
                    args...) where {M<:AbstractModel,S, T<:Real}

    @inbounds for i in eachindex(u,v,vis)
        vis[i] = visibility(m, u[i], v[i], args...)
    end
    return vis
end

function visibilities(m::M,
                      u::AbstractVector{T},
                      v::AbstractVector{T},
                      args...) where {S,M<:AbstractModel{S}, T<:Real}
    vis = StructArray{Complex{S}}(re=similar(u,S), im=similar(v,S))
    return visibilities!(vis, m, u, v, args...)
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
    intensity(model::AbstractModel, x, y, args...)
Returns the intensity for the model at x and y with args...
"""
function intensity end

include("stokesmatrix.jl")
include("geometric_models.jl")
include("radio_image_models.jl")
include("modifiers.jl")
