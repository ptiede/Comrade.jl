export DeltaPulse, SqExpPulse, BSplinePulse

"""
    $(TYPEDEF)
A dirac comb pulse function. This means the image is just the
dicrete Fourier transform
"""
struct DeltaPulse{T} <: Pulse end
DeltaPulse() = DeltaPulse{Float64}()
# This should really be a delta function but whatever
@inline κ(::DeltaPulse{T}, x) where {T} = one(T)
@inline ω(::DeltaPulse{T}, u) where {T} = one(T)


"""
    $(TYPEDEF)
Normalized square exponential kernel, i.e. a Gaussian. Note the
smoothness is modfied with `ϵ` which is the inverse variance in units of
1/pixels².
"""
struct SqExpPulse{T} <: Pulse
    ϵ::T
end
@inline @fastmath κ(b::SqExpPulse, x) = exp(-0.5*b.ϵ^2*x^2)/sqrt(2*π/b.ϵ^2)
@inline κflux(::SqExpPulse{T}) where {T} = one(T)
@inline @fastmath ω(b::SqExpPulse, u) = exp(-2*(π*u/b.ϵ)^2)


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
@inline ω(::BSplinePulse{N}, u) where {N} = sinc(π*u)^N
@inline κflux(::BSplinePulse) = 1.0

@inline κ(b::BSplinePulse{0}, x::T) where {T} = abs(x) < 0.5 ? one(T) : zero(T)

@inline function κ(b::BSplinePulse{1}, x::T) where {T}
    mag = abs(x)
    return mag < 1 ? 1-mag : zero(T)
end

@inline function κ(b::BSplinePulse{3}, x::T) where {T}
    mag = abs(x)
    if mag < 1
        return evalpoly(mag, (4, 0, -6, 3))/6
    elseif 1 ≤ mag < 2
        return evalpoly(mag, (8, -12, 6, -1))/6
    else
        return zero(T)
    end
end
