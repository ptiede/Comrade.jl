using MeasureBase
using SpecialFunctions: besseli
using Random
using KeywordCalls
import Statistics
import Distributions as Dists

export Rice, CPVonMises, CPNormal, CMvNormal, AmpNormal

@parameterized Rice(ν, σ) ≃ Lebesgue(ℝ₊)

@kwstruct Rice(ν, σ)

@parameterized ComplexNormal(μ, σ)
@kwstruct ComplexNormal(μ, σ)
@kwstruct ComplexNormal(μ, τ)

function MeasureBase.logdensity(d::ComplexNormal{(:μ, :σ)}, x)
    sum = zero(eltype(d.σ))
    @inbounds @fastmath for i in eachindex(x)
        sum += abs2((x[i] - d.μ[i])/d.σ[i])
    end
    #return -sum(abs2, (x .- d.μ)./d.σ)/2
    return -sum/2
end

function MeasureBase.logdensity(d::ComplexNormal{(:μ, :τ)}, x)
    sum = zero(eltype(real(d.μ)))
    @inbounds @fastmath for i in eachindex(x)
        sum += abs2((x[i] - d.μ[i])*d.τ[i])
    end
    #return -sum(abs2, (x .- d.μ).*d.τ)/2
    return -sum/2
end

function Base.rand(rng::AbstractRNG, T::Type, d::ComplexNormal{(:μ, :σ)})
    x1 = randn(rng, T)*d.σ + real(d.μ)
    x2 = randn(rng, T)*d.σ + imag(d.μ)
    return x1 + 1im*x2
end

function Base.rand(rng::AbstractRNG, T::Type, d::ComplexNormal{(:μ, :τ)})
    x1 = randn(rng, T)/d.τ + real(d.μ)
    x2 = randn(rng, T)/d.τ + imag(d.μ)
    return x1 + 1im*x2
end

function MeasureBase.logdensity(d::Rice{(:ν, :σ)}, x)
    li0 = log(besselix(0, x*d.ν/d.σ^2)) + abs(x*d.ν/d.σ^2)
    log(x/d.σ^2) + li0 - (x^2 + d.ν^2)/(2*d.σ^2)
end

function Base.rand(rng::AbstractRNG, T::Type, d::Rice{(:ν, :σ)})
    x1 = randn(rng, T)*d.σ + d.ν
    x2 = randn(rng, T)*d.σ
    return hypot(x1,x2)
end

#### Statistics
L12(x) = exp(x/2)*( (1-x)*besseli(0,-x/2) - x*besseli(1,-x/2) )
Statistics.mean(d::Rice{(:ν,:σ)}) = d.σ*sqrt(π/2)*L12(-0.5*(d.ν/(d.σ))^2)
Statistics.var(d::Rice{(:ν, :σ)}) = 2*d.σ^2 + d.ν^2 - π*d.σ^2/2*L12(-0.5*(d.ν/(d.σ))^2)^2

@parameterized CPVonMises(μ, κ)

@kwstruct CPVonMises(μ, κ)
@kwstruct CPVonMises(μ, σ)
const log2π = log(2π)
function MeasureBase.logdensity(d::CPVonMises{(:μ, :κ)}, x)
    T = eltype(d.μ)
    sum = zero(T)
    @inbounds for i = eachindex(x)
        sum += d.κ[i]*(cos(x[i]-d.μ[i])-1) #- log(besselix(0.0, d.κ[i])) - log2π
    end
    return sum
end

function MeasureBase.logdensity(d::CPVonMises{(:μ, :σ)},x)
    sum = zero(eltype(d.μ))
    @inbounds for i = eachindex(x)
        sum += (cos(x[i]-d.μ[i])-1)/d.σ[i]^2 #- log(besselix(0.0, d.κ[i])) - log2π
    end
    return sum
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPVonMises{(:μ, :κ)})
    d = Dists.VonMises.(d.μ, d.κ)
    return rand.(Ref(rng), Ref(T), d)
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPVonMises{(:μ, :σ)})
    d = @. Dists.VonMises(d.μ, 1/d.σ^2)
    return rand.(Ref(rng), Ref(T), d)
end

@parameterized AmpNormal(μ, τ)

@kwstruct AmpNormal(μ, τ)

function MeasureBase.logdensity(d::AmpNormal{(:μ, :τ)}, x)
    T = eltype(d.μ)
    sum = zero(T)
    @inbounds for i = eachindex(x)
        sum += -(d.τ[i]*(d.μ[i] - x[i]))^2/2.0
    end
    return sum
end

function Base.rand(rng::AbstractRNG, T::Type, d::AmpNormal{(:μ, :τ)})
    nd = length(d.μ)
    return randn(rng, T, nd)./d.τ + d.μ
end

@parameterized CPNormal(μ, σ)

Base.minimum(::CPNormal) = -Inf
Base.maximum(::CPNormal) = Inf

const log2π = log(2π)
function MeasureBase.logdensity(dist::CPNormal{(:μ, :σ)}, x::Real)
    μ,σ = dist.μ, dist.σ
    s,c = sincos(x)
    dθ = atan(s, c)
    return  -(abs2(dθ/σ) + log2π)/2
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPNormal)
    return d.μ + d.σ*randn(rng, T)
end
