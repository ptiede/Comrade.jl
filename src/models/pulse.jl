export DeltaPulse, SqExpPulse, BSplinePulse, BicubicPulse, RaisedCosinePulse

"""
    Pulse
Pixel response function for a radio image model. This makes
a discrete sampling continuous by picking a certain *smoothing*
kernel for the image.

# Notes
To see the implemented Pulses please use the subtypes function i.e.
`subtypes(Pulse)`
"""
abstract type Pulse <: AbstractModel end

visanalytic(::Type{<:Pulse}) = IsAnalytic()
imanalytic(::Type{<:Pulse}) = IsAnalytic()
isprimitive(::Type{<:Pulse}) = IsPrimitive()

@inline intensity_point(p::Pulse, ps) = κ(p::Pulse, ps.X)*κ(p::Pulse, ps.Y)
@inline visibility_point(p::Pulse, u, v, time, freq) = ω(p::Pulse, u)*ω(p::Pulse, v)


flux(p::Pulse) = κflux(p)^2

"""
    $(TYPEDEF)
A dirac comb pulse function. This means the image is just the
dicrete Fourier transform
"""
struct DeltaPulse{T} <: Pulse end
DeltaPulse() = DeltaPulse{Float64}()
# This should really be a delta function but whatever
κflux(::DeltaPulse{T}) where {T} = one(T)
@inline κ(::DeltaPulse{T}, x) where {T} = abs(x) < 0.5 ? one(T) : zero(T)
@inline ω(::DeltaPulse{T}, u) where {T} = complex(one(T))
@inline radialextent(::Pulse) = 1.0


@doc raw"""
    $(TYPEDEF)
Uses the basis spline (BSpline) kernel of order `N`. These are the kernel that come
from recursively convolving the tophat kernel
```math
    B_0(x) = \begin{cases} 1 & |x| < 1 \\ 0 & otherwise \end{cases}
```
`N` times.

## Notes

BSpline kernels have a number of nice properties:
1. Simple frequency response $\sinc(u/2)^N$
2. preserve total intensity

For `N`>1 these kernels aren't actually interpolation kernels however, this doesn't matter
for us.

Currently only the 0,1,3 order kernels are implemented.
"""
struct BSplinePulse{N} <: Pulse end
@inline ω(::BSplinePulse{N}, u) where {N} = complex(sinc(u)^(N+1))
@inline κflux(::BSplinePulse) = 1.0

@inline κ(::BSplinePulse{0}, x::T) where {T} = abs(x) < 0.5 ? one(T) : zero(T)

@inline function κ(::BSplinePulse{1}, x::T) where {T}
    mag = abs(x)
    return mag < 1 ? 1-mag : zero(T)
end

@inline function κ(::BSplinePulse{3}, x::T) where {T}
    mag = abs(x)
    if mag < 1
        return evalpoly(mag, (4, 0, -6, 3))/6
    elseif 1 ≤ mag < 2
        return evalpoly(mag, (8, -12, 6, -1))/6
    else
        return zero(T)
    end
end

"""
    $(TYPEDEF)
    BicubicPulse(b = 0.5)

The bicubic pulse for imaging. This pulse tends to have a flat spectrum but for most values
of `b` can produce negative intensities in an image. This is the pulse used in
[`Broderick et al. 2020`](https://iopscience.iop.org/article/10.3847/1538-4357/ab9c1f).
"""
struct BicubicPulse{T} <: Pulse
    b::T
end

BicubicPulse() = BicubicPulse{Float64}(-0.5)

function κ(k::BicubicPulse, x::T) where {T}
    mag = abs(x)
    b = k.b
    if mag < 1
        return evalpoly(mag, (one(T), zero(T), -(b+3), b+2))
    elseif 1 ≤ mag < 2
        return b*evalpoly(mag, (-T(4), T(8), -T(5), one(T)))
    else
        return zero(T)
    end
end

function ω(m::BicubicPulse, u)
    b = m.b
    k = 2π*u
    abs(k) < 1e-2 && return 1 - (2*b + 1)*k^2/15 + (16*b + 1)*k^4/560 + 0im
    s,c = sincos(k)
    k3 = k^3
    k4 = k3*k
    c2 = c^2 - s^2
    return -4*s*(2*b*c + 4*b + 3)*inv(k3) + 12*inv(k4)*(b*(1-c2) + 2*(1-c)) + 0im
end

"""
    $(TYPEDEF)
    RaisedCosinePulse()
    RaisedCosinePulse(rolloff)

Raised cosine pulse function. This tends to be a very flat response, where the roll off
controls the speed of decay. By default we set `rolloff = 0.5`.
"""
struct RaisedCosinePulse{T} <: Pulse
    rolloff::T
end

κflux(::Pulse) = 1

RaisedCosinePulse() = RaisedCosinePulse{Float64}(0.5)

function κ(k::RaisedCosinePulse, x::T) where {T}
    mag = abs(x)
    β = k.rolloff
    if 2*mag < 1-β
        return one(T)
    elseif 1-β <= 2*mag <= 1+β
        return 1/2*(1 + cospi((mag - (1-β)/2)/β))
    else
        return zero(T)
    end
end

function ω(k::RaisedCosinePulse, u::T) where {T}
    β = k.rolloff
    abs(u) ≈ inv(2*β) && return complex(π/4*sinc(inv(2*β)))
    return complex(sinc(u)*cospi(β*u)*inv(1 - (2β*u)^2))
end
