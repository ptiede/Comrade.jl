import Distributions
using Random
using SpecialFunctions: besseli

struct CPNormal{T,S} <: Distributions.ContinuousUnivariateDistribution
    μ::T
    σ::S
end

Base.minimum(::CPNormal) = -Inf
Base.maximum(::CPNormal) = Inf

const log2π = log(2*π)
function Distributions.logpdf(dist::CPNormal, x::Real)
    μ,σ = dist.μ, dist.σ
    s,c = sincos(x-μ)
    dθ = atan(s, c)
    return  -(abs2(dθ/σ) + log2π)/2 - log(σ)
end

function Base.rand(rng::AbstractRNG, d::CPNormal)
    return d.μ + d.σ*randn(rng)
end

Base.rand(rng::AbstractRNG, ::Type{Float64}, d::CPNormal) = rand(rng, d)


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
