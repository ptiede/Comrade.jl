using TupleVectors
using AbstractMCMC

using AbstractDifferentiation

struct IsCube end

struct IsFlat end

struct HasDeriv end
struct NoDeriv end

export sample


samplertype(::Type) = ArgumentError("samplertype not specified")

include(joinpath(@__DIR__, "pullbacks.jl"))

function AbstractMCMC.sample(post::Posterior, sampler::S, args...; init_params=nothing, kwargs...) where {S}
    θ0 = init_params
    if isnothing(init_params)
        θ0 = prior_sample(post)
    end
    return _sample(samplertype(S), post, sampler, args...; init_params=θ0, kwargs...)
end

function _sample(::IsFlat, post, sampler, args...; init_params, kwargs...)
    tpost = asflat(post)
    θ0 = HypercubeTransform.inverse(tpost, init_params)
    return sample(tpost, sampler, args...; init_params=θ0, kwargs...)
end

function _sample(::IsCube, post, sampler, args...; init_params, kwargs...)
    tpost = ascube(post)
    θ0 = HypercubeTransform.inverse(tpost, init_params)
    return sample(tpost, sampler, args...; init_params=θ0, kwargs...)
end
