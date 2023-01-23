export Gaussian,
        Disk,
        MRing,
        Crescent,
        ConcordanceCrescent,
        ExtendedRing,
        Ring,
        ParabolicSegment,
        Wisp,
        Butterworth,
        SlashedDisk

# helper functions for below
@inline _getuv(p) = (p.U, p.V)
@inline _getxy(p) = (p.X, p.Y)


"""
$(TYPEDEF)
A type that defines it is a geometric model. These are usually
primitive models, and are usually analytic in Fourier and the image domain.
As a result a user only needs to implement the following methods

- `visibility_point`
- `intensity_point`
- `radialextent`

Note that if the geometric model isn't **analytic** then the usual methods listed
in [`Comrade.AbstractModel`](@ref) for non-analytic models need to be implemented.
"""
abstract type GeometricModel <: AbstractModel end
@inline flux(::GeometricModel) = 1.0

@inline isprimitive(::Type{<:GeometricModel}) = IsPrimitive()

@inline visanalytic(::Type{<:GeometricModel}) = IsAnalytic()
@inline imanalytic(::Type{<:GeometricModel}) = IsAnalytic()



"""
    $(TYPEDEF)

Gaussian with unit standard deviation and flux.

By default if T isn't given, `Gaussian` defaults to `Float64`
"""
struct Gaussian{T} <: GeometricModel end
Gaussian() = Gaussian{Float64}()
radialextent(::Gaussian) = 5.0


@inline function intensity_point(::Gaussian, p)
    x,y = _getxy(p)
    return exp(-(x^2+y^2)/2)/2π
end

@inline function visibility_point(::Gaussian{T}, u, v, time, freq) where {T}
    return exp(-2π^2*(u^2 + v^2)) + zero(T)im
end




@doc raw"""
    Disk{T}() where {T}

Tophat disk geometrical model, i.e. the intensity profile
```math
    I(x,y) = \begin{cases} \pi^{-1} & x^2+y^2 < 1 \\ 0 & x^2+y^2 \geq 0 \end{cases}
```
i.e. a unit radius and unit flux disk.

By default if T isn't given, `Disk` defaults to `Float64`
"""
struct Disk{T} <: GeometricModel end
Disk() = Disk{Float64}()

@inline function intensity_point(::Disk{T}, p) where {T}
    x,y = _getxy(p)
    r = x^2 + y^2
    return r < 1 ?  one(T)/(π) : zero(T)
end

@inline function visibility_point(::Disk{T}, u, v, time, freq) where {T}
    ur = 2π*(hypot(u,v) + eps(T))
    return 2*besselj1(ur)/(ur) + zero(T)im
end

radialextent(::Disk) = 3.0


@doc raw"""
    SlashedDisk{T}(slash::T) where {T}

Tophat disk geometrical model, i.e. the intensity profile
```math
    I(x,y) = \begin{cases} \pi^{-1} & x^2+y^2 < 1 \\ 0 & x^2+y^2 \geq 0 \end{cases}
```
i.e. a unit radius and unit flux disk.

By default if T isn't given, `Disk` defaults to `Float64`
"""
struct SlashedDisk{T} <: GeometricModel
    slash::T
end


function intensity_point(m::SlashedDisk{T}, p) where {T}
    x,y = _getxy(p)
    r2 = x^2 + y ^2
    s = 1 - m.slash
    norm = 2*inv(π)/(1+s)
    if  r2 < 1
        return norm/2*((1+y) + s*(1-y))
    else
        return zero(T)
    end
end

function visibility_point(m::SlashedDisk{T}, u, v, time, freq) where {T}
    k = 2π*sqrt(u^2 + v^2) + eps(T)
    s = 1-m.slash
    norm = 2/(1+s)/k

    b0outer = besselj0(k)
    b1outer = besselj1(k)
    b2outer = besselj(2,k)

    v1 = (1+s)*b1outer
    v3 = -2im*π*u*(1-s)*(b0outer-b2outer-2*b1outer/k)/(2*k)
    return norm*(v1+v3)
end

radialextent(::SlashedDisk) = 3.0


"""
    $(TYPEDEF)

A infinitely thin ring model, whose expression in the image domain is
    I(r,θ) = δ(r - 1)/2π
i.e. a unit radius and flux delta ring.

By default if `T` isn't given, `Gaussian` defaults to `Float64`
"""
struct Ring{T} <: GeometricModel end
Ring() = Ring{Float64}()
radialextent(::Ring) = 3/2

@inline function intensity_point(::Ring{T}, p) where {T}
    x,y = _getxy(p)
    r = hypot(x,y)
    θ = atan(x,y)
    dr = 1e-2
    if (abs(r-1) < dr/2)
        acc = one(T)
        return acc/(2π*dr)
    else
        return zero(T)
    end
end



@inline function visibility_point(::Ring{T}, u, v, time, freq) where {T}
    k = 2π*sqrt(u^2 + v^2) + eps(T)
    vis = besselj0(k) + zero(T)*im
    return vis
end

struct Butterworth{N, T} <: GeometricModel end

"""
    Butterworth{N}()
    Butterworth{N, T}()

Construct a model that corresponds to the Butterworth filter of order `N`.
The type of the output is given by `T` and if not given defaults to `Float64`
"""
Butterworth{N}() where {N} = Butterworth{N,Float64}()

radialextent(b::Butterworth) = 5
flux(::Butterworth{N,T}) where {N,T} = one(T)

visanalytic(::Type{<:Butterworth}) = IsAnalytic()
imanalytic(::Type{<:Butterworth}) = NotAnalytic()

function visibility_point(::Butterworth{N,T}, u, v, time, freq) where {N,T}
    b = hypot(u,v) + eps(T)
    return complex(inv(sqrt(1 + b^(2*N))))
end





"""
    $(TYPEDEF)
m-ring geometric model. This is a infinitely thin unit flux delta ring
whose angular structure is given by a Fourier expansion. That is,

    I(r,θ) = (2π)⁻¹δ(r-1)∑ₙ(αₙcos(nθ) - βₙsin(nθ))

The `N` in the type defines the order of the Fourier expansion.



# Fields
$(FIELDS)
"""
struct MRing{T, V<:Union{AbstractVector{T}, NTuple}} <: GeometricModel
    """
    Real Fourier mode coefficients
    """
    α::V
    """
    Imaginary Fourier mode coefficients
    """
    β::V
    function MRing(α::V, β::V) where {V <: Union{AbstractVector, NTuple}}
        @argcheck length(α) == length(β) "Lengths of real/imag components must be equal in MRing"
        return new{eltype(α), V}(α, β)
    end
end


"""
    MRing(c::AbstractVector{<:Complex})

Construct an MRing geometric model from a complex vector `c`
that correspond to the real and imaginary (or cos and sin) coefficients
of the Fourier expansion. The `N` in the type defines the order of
the Fourier expansion.

"""
function MRing(c::Union{AbstractVector{<:Complex}, NTuple{N,<:Complex}}) where {N}
    α = real.(c)
    β = imag.(c)
    return MRing(α, β)
end

function MRing(a::Real, b::Real)
    aT,bT = promote(a,b)
    return MRing((aT,), (bT,))
end

# Depreciate this method since we are moving to vectors for simplificty
#@deprecate MRing(a::Tuple, b::Tuple) MRing(a::AbstractVector, b::AbstractVector)

radialextent(::MRing) = 1.5


@inline function intensity_point(m::MRing{T}, p) where {T}
    x,y = _getxy(p)
    r = hypot(x,y)
    θ = atan(x,y)
    dr = 0.025
    if (abs(r-1) < dr/2)
        acc = one(T)
        for n in eachindex(m.α, m.β)
            s,c = sincos(n*θ)
            acc += m.α[n]*c - m.β[n]*s
        end
        return acc/(2π*dr)
    else
        return zero(T)
    end
end


@inline function visibility_point(m::MRing{T}, u, v, time, freq) where {T}
    return _mring_vis(m, u, v)
end

@inline function _mring_vis(m::MRing{T}, u, v) where {T}
    (;α, β) = m
    k = 2π*sqrt(u^2 + v^2) + eps(T)
    vis = besselj0(k) + zero(T)*im
    θ = atan(u, v)
    @inbounds for n in eachindex(α, β)
        s,c = sincos(n*θ)
        vis += 2*(α[n]*c - β[n]*s)*(1im)^n*besselj(n, k)
    end
    return vis
end

function _mring_adjoint(α, β, u, v)
    T = eltype(α)
    ρ = hypot(u,v)
    k = 2π*ρ + eps(T)
    θ = atan(u,v)
    vis = complex(besselj0(k))

    j0 = besselj0(k)
    j1 = besselj1(k)
    #bj = Base.Fix2(besselj, k)
    #jn = ntuple(bj, length(α))

    ∂u = -complex(j1*2π*u/ρ)
    ∂v = -complex(j1*2π*v/ρ)
    ∂α = similar(α ,Complex{T})
    ∂β = similar(α ,Complex{T})

    ∂ku = 2π*u/ρ
    ∂kv = 2π*v/ρ
    ∂θu = v/ρ^2
    ∂θv = -u/ρ^2

    for n in eachindex(α, β)
        s,c = sincos(n*θ)
        imn = (1im)^n
        jn = besselj(n,k)
        ∂α[n] =  2*c*jn*imn
        ∂β[n] = -2*s*jn*imn
        dJ = j0 - n/k*jn
        j0 = jn

        visargc = 2*imn*(-α[n]*s - β[n]*c)
        visarg  = 2*imn*(α[n]*c - β[n]*s)
        vis += visarg*jn

        ∂u +=  n*visargc*jn*∂θu + visarg*dJ*∂ku
        ∂v +=  n*visargc*jn*∂θv + visarg*dJ*∂kv

    end
    return vis, ∂α, ∂β, ∂u, ∂v
end

function ChainRulesCore.rrule(::typeof(_mring_vis), m::MRing, u, v)
    (;α, β) = m
    vis, ∂α, ∂β, ∂u, ∂v = _mring_adjoint(α, β, u, v)

    function _mring_pullback(Δv)
        return (NoTangent(), Tangent{typeof(m)}(α=real(Δv'.*∂α), β=real(Δv'.*∂β)), real(Δv'.*∂u), real(Δv'.*∂v))
    end

    return vis, _mring_pullback

end

"""
    $(TYPEDEF)

Creates a [Kamruddin and Dexter](https://academic.oup.com/mnras/article/434/1/765/1005984)
crescent model. This works by composing two disk models together.

# Arguments
- `router`: The radius of the outer disk
- `rinner`: The radius of the inner disk
- `shift`: How much the inner disk radius is shifted (positive is to the right)
- `floor`: The floor of the inner disk 0 means the inner intensity is zero and 1 means it is a large disk.
"""
function Crescent(router, rinner, shift, floor)
    m = stretched(Disk(), router, router)*(π*router^2) - shifted(stretched(Disk(), rinner, rinner)*((1-floor)*π*rinner^2), shift, zero(typeof(shift)))
    return m/flux(m)
end

"""
    $(TYPEDEF)
Creates the ConcordanceCrescent model, i.e. a flat-top crescent
with a displacment and a slash and shadow depth.
Note this creates a crescent with unit flux.
If you want a different flux please use the `renomed`
modifier.

## Fields
$(FIELDS)

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

function intensity_point(m::ConcordanceCrescent{T}, p) where {T}
    x,y = _getxy(p)
    r2 = x^2 + y ^2
    norm = _crescentnorm(m)
    if (r2 < m.router^2 && (x-m.shift)^2 + y^2 > m.rinner^2 )
        return norm/2*((1+x/m.router) + m.slash*(1-x/m.router))
    else
        return zero(T)
    end
end

function visibility_point(m::ConcordanceCrescent{T}, u, v, time, freq) where {T}
    k = 2π*sqrt(u^2 + v^2) + eps(T)
    norm = π*_crescentnorm(m)/k
    phaseshift = cispi(2*m.shift*u)
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

The formula in the image domain is given by

    I(r,θ) = βᵅrᵅ⁻²exp(-β/r)/2πΓ(α)

where `α = shape` and `β = shape+1`

# Note
We mainly use this as an example of a non-analytic Fourier transform
(although it has a complicated expression)

# Fields
$(FIELDS)

Note that if `T` isn't specified at construction then it defaults to `Float64`.
"""
struct ExtendedRing{F<:Number} <: GeometricModel
    """shape of the radial distribution"""
    shape::F
end
visanalytic(::Type{<:ExtendedRing}) = NotAnalytic()

radialextent(m::ExtendedRing) = 6.0

function intensity_point(m::ExtendedRing, p)
    x,y = _getxy(p)
    r = hypot(x, y) + eps()
    β = (m.shape + 1)
    α = m.shape
    β^α*r^(-α-2)*exp(-β/r)/gamma(α)/(2*π)
end


"""
    $(TYPEDEF)

A infinitely thin parabolic segment in the image domain.
The segment is centered at zero, with roots ±1 and a yintercept of 1.

Note that if `T` isn't specified at construction then it defaults to `Float64`.
"""
struct ParabolicSegment{T} <: GeometricModel end
ParabolicSegment() = ParabolicSegment{Float64}()
radialextent(::ParabolicSegment{T}) where {T} = one(T)*sqrt(2)

"""
    ParabolicSegment(a::Number, h::Number)

A parabolic segment with x-intercepts `±a` and a yintercept of `h`.

# Note
This is just a convenience function for `stretched(ParabolicSegment(), a, h)`
"""
@inline function ParabolicSegment(a::Number, h::Number)
    # Define stretched model from unital model
    stretched(ParabolicSegment(), a, h)
end

function intensity_point(::ParabolicSegment{T}, p) where {T}
    x,y = _getxy(p)
    yw = (1-x^2)
    if abs(y - yw) < 0.01/2 && abs(x) < 1
        return 1/(2*0.01)
    else
        return zero(T)
    end
end

function visibility_point(::ParabolicSegment{T}, u, v, time, freq) where {T}
    ϵ = sqrt(eps(T))
    vϵ = v + ϵ + 0im
    phase = cispi(3/4 + 2*vϵ + u^2/(2vϵ))
    Δ1 = erf(√(π/(2vϵ))*cispi(1/4)*(u-2vϵ))
    Δ2 = erf(√(π/(2vϵ))*cispi(1/4)*(u+2vϵ))
    return phase/(√(2vϵ))*(Δ1-Δ2)/4
end
