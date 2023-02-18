module ComradeDynesty

using Comrade

using AbstractMCMC
using InferenceObjects
using Reexport
using Random
using StatsBase
using Tables

@reexport using Dynesty


Comrade.samplertype(::Type{<:NestedSampler}) = Comrade.IsCube()
Comrade.samplertype(::Type{<:DynamicNestedSampler}) = Comrade.IsCube()

"""
    AbstractMCMC.sample(post::Comrade.Posterior, smplr::Dynesty.NestedSampler, args...; kwargs...)
    AbstractMCMC.sample(post::Comrade.Posterior, smplr::Dynesty.DynamicNestedSampler, args...; kwargs...)

Sample the posterior `post` using `Dynesty.jl` `NestedSampler/DynamicNestedSampler` sampler.
The `args/kwargs`
are forwarded to `Dynesty` for more information see its [docs](https://github.com/ptiede/Dynesty.jl)

This returns a tuple where the first element is a InferenceData object that contains the weighted
samples produced by dynesty.

To create equally weighted samples the user can do
```julia
samples, dy = sample(post, NestedSampler(dimension(post), 1000))
equal_weighted_chain = ComradeDynesty.equalresample(samples, 1000)
```
"""
function AbstractMCMC.sample(::Random.AbstractRNG, post::Comrade.TransformedPosterior,
                             sampler::Union{NestedSampler, DynamicNestedSampler}
                             ; init_params=nothing,
                             kwargs...)
    ℓ = logdensityof(post)
    kw = delete!(Dict(kwargs), :init_params)
    res = dysample(ℓ, identity, sampler; kw...)
    # Make sure that res["sample"] is an array and use transpose
    samples, weights = transpose(Dynesty.PyCall.PyArray(res["samples"])), exp.(res["logwt"] .- res["logz"][end])
    chain = [transform.(Ref(post), eachcol(samples))]
    sample_stats = (log_density = [Array(Dynesty.PyCall.PyArray(res["logl"]))],
                     weights = [weights],
                    )
    metadata = Dict(:logz => res["logz"][end], :logzerr => res["logzerr"][end])
    library = "Comrade Dynesty.jl"
    return from_namedtuple(chain; sample_stats, library, attrs=metadata), res
end

function equalresample(res::InferenceData, nsamples)
    weights = vec(res.sample_stats.weights)
    eq = sample(Tables.rowtable(res.posterior[chain=1]), Weights(weights), nsamples)
    return from_namedtuple([eq])
end



end
