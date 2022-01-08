using MeasureBase
using SpecialFunctions: besseli
using Random
using KeywordCalls
import Statistics: mean, var

export Rice, CPVonMises

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
    return -
end

struct CPVonMises{T,S} <: Distributions.ContinuousUnivariateDistribution
    μ::T
    σ::S
    I0κx::S
end

function CPVonMises(μ, σ)
    CPVonMises(μ, σ, besselix(zero(typeof(σ)), 1/σ^2))
end

function Distributions.logpdf(dist::CPVonMises, x::Real)
    μ,σ = dist.μ, dist.σ
    dθ = (cos(x-μ)-1)/σ^2
    return dθ - log(dist.I0κx) - log2π
end

Base.minimum(::CPVonMises) = -Inf
Base.maximum(::CPVonMises) = Inf


function Base.rand(rng::AbstractRNG, d::CPVonMises)
    return d.μ + d.σ*randn(rng)
end

Base.rand(rng::AbstractRNG, ::Type{Float64}, d::CPVonMises) = rand(rng, d)
