using TupleVectors
using AbstractMCMC

using AbstractDifferentiation

struct IsCube end

struct IsFlat end

struct HasDeriv end
struct NoDeriv end

export sample



include(joinpath(@__DIR__, "pullbacks.jl"))

function AbstractMCMC.sample(post::Posterior, sampler::S, args...; init_params=nothing, kwargs...) where {S}
    return _sample(samplertype(S), post, sampler, args...; init_params, kwargs...)
end

function _sample(::IsFlat, post, sampler, args...; init_params, kwargs...)
    tpost = asflat(post)
    return sample(tpost, sampler, args...; init_params, kwargs...)
end

function _sample(::IsCube, post, sampler, args...; init_params, kwargs...)
    tpost = ascube(post)
    return sample(tpost, sampler, args...; init_params, kwargs...)
end
