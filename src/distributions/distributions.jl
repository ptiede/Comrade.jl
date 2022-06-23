using SpecialFunctions: besseli
using Random
using KeywordCalls
import Statistics
import Distributions as Dists

export CPVonMises, CPNormal, CMvNormal, AmpNormal


"""
    ComplexNormal

Uncorrelated complex Normal Measure. This is the default likelihood used for visibilties.
"""
ComplexNormal

MT.@parameterized ComplexNormal(μ, σ)
MT.@kwstruct ComplexNormal(μ, σ)
MT.@kwstruct ComplexNormal(μ, τ)

MB.basemeasure(::ComplexNormal) = MT.Lebesgue()

function MB.logdensity_def(d::ComplexNormal{(:μ, :σ)}, x)
    #sum = zero(eltype(d.σ))
    #@inbounds @fastmath for i in eachindex(x)
    #    sum += abs2((x[i] - d.μ[i])/d.σ[i])
    #end
    return -sum(abs2, (x .- d.μ)./d.σ)/2
    #return -sum/2
end

MB.insupport(::ComplexNormal, x) = true
MB.insupport(::ComplexNormal) = MB.Returns(true)


function MB.logdensity_def(d::ComplexNormal{(:μ, :τ)}, x)
    #sum = zero(eltype(real(d.μ)))
    #@inbounds @fastmath for i in eachindex(x)
    #    sum += abs2((x[i] - d.μ[i])*d.τ[i])
    #end
    return -sum(abs2, (x .- d.μ).*d.τ)/2
    #return -sum/2
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


"""
    CPVonMises

The von Mises distribution used for closure phases. Typically μ and κ=1/σ² are vectors
to allow for simpler derivative rules.
"""
CPVonMises

MT.@parameterized CPVonMises(μ, κ)

MT.@kwstruct CPVonMises(μ, κ)
MT.@kwstruct CPVonMises(μ, σ)

MB.basemeasure(::CPVonMises) = MT.Lebesgue()

const log2π = log(2π)
function MB.logdensity_def(d::CPVonMises{(:μ, :κ)}, x)
    #T = eltype(d.μ)
    #sum = zero(T)
    #@inbounds for i = eachindex(x)
    #    sum += d.κ[i]*(cos(x[i]-d.μ[i])-1) #- log(besselix(0.0, d.κ[i])) - log2π
    #end
    #return sum
    dθ = @. d.κ*(cos(x - d.μ) - 1)
    return sum(dθ)
end

MB.insupport(::CPVonMises, x) = true
MB.insupport(::CPVonMises) = MB.Returns(true)


function MB.logdensity_def(d::CPVonMises{(:μ, :σ)},x)
    #sum = zero(eltype(d.μ))
    #@inbounds for i = eachindex(x)
    #    sum += (cos(x[i]-d.μ[i])-1)/d.σ[i]^2 #- log(besselix(0.0, 1/d.σ[i]^2)) - log2π
    #end
    dθ = @. inv(d.σ)^2*(cos(x-d.μ) - 1)
    return sum(dθ)
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPVonMises{(:μ, :κ)})
    d = Dists.VonMises.(d.μ, d.κ)
    return rand.(Ref(rng), Ref(T), d)
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPVonMises{(:μ, :σ)})
    d = @. Dists.VonMises(d.μ, 1/d.σ^2)
    return rand.(Ref(rng), Ref(T), d)
end


"""
    AmpNormal

Visibility amplitude likelihood distribution. Typically μ, τ=1/σ are vectors to allow for
more efficient derivatives.
"""
AmpNormal


MT.@parameterized AmpNormal(μ, τ)
MT.@kwstruct AmpNormal(μ, τ)
MB.basemeasure(::AmpNormal) = MT.Lebesgue()


function MB.logdensity_def(d::AmpNormal{(:μ, :τ)}, x)
    #T = eltype(d.μ)
    #sum = zero(T)
    #@inbounds for i = eachindex(x)
    #    sum += -(d.τ[i]*(d.μ[i] - x[i]))^2/2.0 #-0.5*log2π + log(d.τ[i])
    #end
    return -sum(abs2, d.τ.*(d.μ .- x))/2
end

function Base.rand(rng::AbstractRNG, T::Type, d::AmpNormal{(:μ, :τ)})
    nd = length(d.μ)
    return randn(rng, T, nd)./d.τ + d.μ
end

MB.insupport(::AmpNormal, x) = true
MB.insupport(::AmpNormal) = MB.Returns(true)


MT.@parameterized CPNormal(μ, σ)

Base.minimum(::CPNormal) = -Inf
Base.maximum(::CPNormal) = Inf
MB.basemeasure(::CPNormal) = MT.Lebesgue()

function MB.logdensity_def(dist::CPNormal{(:μ, :σ)}, x::Real)
    μ,σ = dist.μ, dist.σ
    s,c = sincos(x-μ)
    dθ = atan(s, c)
    return  -(abs2(dθ/σ) + log2π)/2
end

function Base.rand(rng::AbstractRNG, T::Type, d::CPNormal)
    return d.μ + d.σ*randn(rng, T)
end
