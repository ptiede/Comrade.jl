export Gaussian, Disk, MRing, Crescent, ConcordanceCrescent, ExtendedRing, Ring

"""
$(TYPEDEF)
A type that defines it is a geometric model. These are usually
primitive models, and are usually analytic in Fourier and the image domain.
"""
abstract type GeometricModel <: AbstractModel end
@inline flux(::GeometricModel) = 1.0

@inline isprimitive(::Type{<:GeometricModel}) = IsPrimitive()

@inline visanalytic(::Type{<:GeometricModel}) = IsAnalytic()
@inline imanalytic(::Type{<:GeometricModel}) = IsAnalytic()



"""
    $(TYPEDEF)
Gaussian geometrical model.
This is a Gaussian with unit flux and standard deviation.

## Notes
To change the Gaussian flux, and shape please use the modifier functions
"""
struct Gaussian{T} <: GeometricModel end
Gaussian() = Gaussian{Float64}()
radialextent(::Gaussian) = 5.0


@inline function intensity_point(::Gaussian, x,y)
    return exp(-(x^2+y^2)/2)/2π
end

@inline function visibility_point(::Gaussian{T}, u, v, args...) where {T}
    return exp(-2π^2*(u^2 + v^2)) + zero(T)im
end




raw"""
    $(TYPEDEF)
Tophat disk geometrical model. The model is given by
```math
    I(x,y) = \begin{cases} \pi^{-1} & x^2+y^2 < 1 \\ 0 & x^2+y^2 \geq 0 \end{cases}
```
"""
struct Disk{T} <: GeometricModel end
Disk() = Disk{Float64}()

@inline function intensity_point(::Disk{T}, x, y, args...) where {T}
    r = x^2 + y^2
    return r < 1 ?  one(T)/(π) : zero(T)
end

@inline function visibility_point(::Disk{T}, x, y, args...) where {T}
    ur = 2π*(hypot(x,y) + eps(T))
    return 2*besselj1(ur)/(ur) + zero(T)im
end

radialextent(::Disk) = 3.0

"""
    $(TYPEDEF)
m-ring geometric model. This corresponds to a delta ring with a fourier expansion
in θ. The m in m-ring refers to the order of the Fourier expansion. The radius is unity.
"""
struct Ring{T} <: GeometricModel end
Ring() = Ring{Float64}()
radialextent(::Ring) = 1.5

@inline function intensity_point(m::Ring, x::Number, y::Number)
    r = hypot(x,y)
    θ = atan(x,y)
    dr = 0.1
    if (abs(r-1) < dr)
        acc = one(T)
        return acc/(2π*dr)
    else
        return zero(T)
    end
end



@inline function visibility_point(m::Ring, u, v, args...)
    k = 2π*sqrt(u^2 + v^2) + eps()
    vis = besselj0(k) + zero(typeof(u))*im
    return vis
end



"""
    $(TYPEDEF)
m-ring geometric model. This corresponds to a delta ring with a fourier expansion
in θ. The m in m-ring refers to the order of the Fourier expansion. The radius is unity.
"""
struct MRing{T,N} <: GeometricModel
    """
    Real Fourier mode coefficients
    """
    α::NTuple{N,T}
    """
    Imaginary Fourier mode coefficients
    """
    β::NTuple{N,T}
end

function MRing{N}(α::T, β::T) where {N,T<:AbstractVector}
    (N != length(α)) && throw("N must be the length of the vector")
    S = promote_type(eltype(α), eltype(β))
    return MRing{S, N}(NTuple{N,S}(α), NTuple{N,S}(β))
end

radialextent(::MRing) = 1.5


@inline function intensity_point(m::MRing{T,N}, x::Number, y::Number) where {T,N}
    r = hypot(x,y)
    θ = atan(x,y)
    dr = 0.1
    if (abs(r-1) < dr)
        acc = one(T)
        for n in 1:N
            s,c = sincos(n*θ)
            acc += m.α[n]*c - m.β[n]*s
        end
        return acc/(2π*dr)
    else
        return zero(T)
    end
end



@inline function visibility_point(m::MRing{T,N}, u, v, args...) where {T,N}
    k = 2π*sqrt(u^2 + v^2) + eps(T)
    vis = besselj0(k) + zero(T)*im
    θ = atan(u, v)
    @inbounds for n in 1:N
        s,c = sincos(n*θ)
        vis += 2*(m.α[n]*c - m.β[n]*s)*(1im)^n*besselj(n, k)
    end
    return vis
end

"""
    $(TYPEDEF)
Creates a [Kamruddin and Dexter](https://academic.oup.com/mnras/article/434/1/765/1005984)
crescent model. This works by composing two disk models together.

# Arguments
- router: The radius of the outer disk
- rinner: The radius of the inner disk
- shift: How much the inner disk radius is shifted (positive is to the right)
- floor: The floor of the inner disk 0 means the inner intensity is zero and 1 means it is a large disk.
"""
function Crescent(router, rinner, shift, floor)
    m = stretched(Disk(), router, router)*(π*router^2) - shifted(stretched(Disk(), rinner, rinner)*((1-floor)*π*rinner^2), shift, zero(typeof(shift)))
    return m/flux(m)
end

"""
    $(TYPEDEF)
Creates the ConcordanceCrescent model, i.e. a flat-top crescent
with a displacment and a slash. Note this creates a crescent with
unit flux. If you want a different flux please use the `renomed`
modifier.

## Fields
## Notes
Unlike the Gaussian and Disk models this does not create the
unit version. In fact, this model could have been created using
the `Disk` and primitives by using Comrade.jl's model composition
functionality.
"""
struct ConcordanceCrescent{T} <: GeometricModel
    """
    Outer radius of the crescent
    """
    router::T
    """
    Inner radius of the crescent
    (i.e. inside this radius there is a hole)
    """
    rinner::T
    """
    Displacment of the inner disk radius
    """
    shift::T
    """
    Strength of the linear slash. Note that
    s∈[0.0,1.0] to ensure positivity in the image.
    """
    slash::T
end

radialextent(m::ConcordanceCrescent) = m.router*1.5

# Crescent normalization to ensure the
function _crescentnorm(m::ConcordanceCrescent)
    f = (1+m.slash)*(m.router^2 - m.rinner^2) -
        (1-m.slash)*m.shift*m.rinner*m.rinner/m.router
    return 2/π/f
end

function intensity_point(m::ConcordanceCrescent{T}, x, y, args...) where {T}
    r2 = x^2 + y ^2
    norm = _crescentnorm(m)
    if (r2 < m.router^2 && (x-m.shift)^2 + y^2 > m.rinner^2 )
        return norm/2*((1+x/m.router) + m.slash*(1-x/m.router))
    else
        return zero(T)
    end
end

function visibility_point(m::ConcordanceCrescent{T}, u, v, args...) where {T}
    k = 2π*sqrt(u^2 + v^2) + eps(T)
    norm = π*_crescentnorm(m)/k
    phaseshift = exp(2im*π*m.shift*u)
    b0outer,b0inner = besselj0(k*m.router), besselj0(k*m.rinner)
    b1outer,b1inner = besselj1(k*m.router), besselj1(k*m.rinner)
    b2outer,b2inner = besselj(2,k*m.router), besselj(2, k*m.rinner)

    v1 = (1+m.slash)*m.router*b1outer
    v2 = ((1+m.slash) + (1-m.slash)*m.shift/m.router)*
            phaseshift*m.rinner*b1inner
    v3 = -2im*π*u*(1-m.slash)*(m.router*b0outer -
                           m.router*b2outer -
                           2*b1outer/k
                          )/(2*k)
    v4 = -2im*π*u*(1-m.slash)*(m.rinner*b0inner -
                          m.rinner*b2inner -
                          2*b1inner/k
                         )/(2*k)*(m.rinner/m.router)*phaseshift
    return norm*(v1-v2+v3-v4)
end

"""
    $(TYPEDEF)
A symmetric extended ring whose radial profile follows an inverse
gamma distributions.

# Note
@e mainly use this as an example of a non-analytic Fourier transform
(although it has a complicated expression)

# Fields

$(FIELDS)

"""
struct ExtendedRing{F} <: GeometricModel
    """radius of peak emission"""
    radius::F
    """shape of the radial distribution"""
    shape::F
end
visanalytic(::Type{<:ExtendedRing}) = NotAnalytic()

radialextent(m::ExtendedRing) = m.radius*5

function intensity_point(m::ExtendedRing, x, y)
    r = hypot(x, y) + eps(m.radius)
    β = m.radius*(m.shape + 1)
    α = m.shape
    β^α*r^(-α-2)*exp(-β/r)/gamma(α)/(2*π)
end
