"""
    $(TYPEDEF)
Gaussian geometrical model.
This is a Gaussian with unit flux and standard deviation.

## Notes
To change the Gaussian flux, and shape please use the modifier functions
"""
struct Gaussian{T} <: GeometricModel{T} end
Gaussian() = Gaussian{Float64}()

function intensity(::Gaussian, x,y, args...)
    return exp(-(x^2+y^2)/2)/sqrt(2π)
end

function visibility(::IsAnalytic, ::Gaussian, u, v, args...)
    return exp(-2π^2*(u^2 + v^2))
end

flux(::GeometricModel{T}) where {T} = one(T)


raw"""
    $(TYPEDEF)
Tophat disk geometrical model. The model is given by
```math
    I(x,y) = \begin{cases} \pi^{-1} & x^2+y^2 < 1 \\ 0 & x^2+y^2 \geq 0 \end{cases}
```
"""
struct Disk{T} <: GeometricModel{T} end
Disk() = Disk{Float64}()

function intensity(::Disk{T}, x, y, args...) where {T}
    r = x^2 + y^2
    return r < 1 ?  π^(-1)*one(T) : zero(T)
end

function visibility(::IsAnalytic, ::Disk{T}, x, y, args...) where {T}
    ur = hypot(x,y) + eps(T)
    return besselj1(2π*ur)/(π*ur)
end
