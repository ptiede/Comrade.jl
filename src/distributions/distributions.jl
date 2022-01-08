using MeasureBase
using SpecialFunctions: besseli
using Random
using KeywordCalls
import Statistics: mean, var
import Distributions as Dists

export Rice, CPVonMises, CPNormal

@parameterized Rice(ν, σ) ≃ Lebesgue(ℝ₊)

@kwstruct Rice(ν, σ)


function MeasureBase.logdensity(d::Rice{(:ν, β)}, x)
    log(x/d.σ^2*besseli(0, x*d.ν/d.σ^2)) - (x^2 + ν^2)/(2*d.σ^2)
end

function Base.rand(rng::AbstractRNG, T::Type, d::Rice{:ν, σ})
    x1 = randn(rng, T)*d.σ + d.ν
    x2 = randn(rng, T)*d.σ
    return hypot(x1,x2)
end

params(d::Rice) = (d.ν, d.σ,)
partype(::Rice{T}) where {T<:Real} = T

#### Statistics
L12(x) = exp(x/2)*( (1-x)*besseli(0,-x/2) - x*besseli(1,-x/2) )
mean(d::Rice{(:ν,:σ)}) = d.σ*Distributions.sqrt2π/2*L12(-(d.ν/(2d.σ))^2)
var(d::Rice{(:ν, :σ)}) = 2*d.σ^2 + d.ν^2 - π*d.σ^2/2*L12(-(ν/(2d.σ))^2)^2

@parameterized CPVonMises(μ, κ)

@kwstruct CPVonMises(μ, κ)

function MeasureBase.logdensity(d::CPVonMises{(:μ, :κ)}, x)
    return d.κ*(cos(x-d.μ)-1)^2
end

function MeasureBase.logdensity(d::CPVonMises{(:μ, :σ)},x)
    return (cos(x-d.μ)-1)^2/d.σ^2
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPVonMises{(:μ, :κ)})
    d = Dists.VonMises(d.μ, d.κ)
    return rand(rng, T, d)
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPVonMises{(:μ, :σ)})
    d = Dists.VonMises(d.μ, 1/d.σ^2)
    return rand(rng, T, d)
end

@parameterized CPNormal(μ, σ)

Base.minimum(::CPNormal) = -Inf
Base.maximum(::CPNormal) = Inf

const log2π = log(2*π)
function MeasureBase.logdensity(dist::CPNormal{(:μ, σ)}, x::Real)
    μ,σ = dist.μ, dist.σ
    s,c = sincos(x)
    dθ = atan(s, c)
    return  -(abs2(dθ/σ) + log2π)/2
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPNormal)
    return d.μ + d.σ*randn(rng, T)
end
