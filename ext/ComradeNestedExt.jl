module ComradeNestedExt

using Comrade
using AbstractMCMC
using NestedSamplers
using Random

"""
    AbstractMCMC.sample(rng::Random.AbstractRNG, post::Comrade.VLBIPosterior, smplr::Nested, args...; kwargs...)

Sample the posterior `post` using `NestedSamplers.jl` `Nested` sampler. The `args/kwargs`
are forwarded to `NestedSampler` for more information see its [docs](https://github.com/TuringLang/NestedSamplers.jl)

This returns a tuple where the first element are the weighted samples from NestedSamplers in a TypedTable.
The second element includes additional information about the samples, like the log-likelihood,
evidence, evidence error, and the sample weights.

To create equally weighted samples the user can use
```julia
using StatsBase
chain, stats = sample(post, NestedSampler(dimension(post), 1000))
equal_weighted_chain = sample(chain, Weights(stats.weights), 10_000)

"""
function AbstractMCMC.sample(rng::Random.AbstractRNG, post::Comrade.VLBIPosterior, sampler::Nested, args...; kwargs...)
    tpost = ascube(post)
    ℓ = logdensityof(tpost)
    model = NestedModel(ℓ, identity)

    samples, stats = sample(rng, model, sampler, args...; chain_type=Array, kwargs...)
    weights = samples[:, end]
    chain = transform.(Ref(tpost), eachrow(samples[:,1:end-1]))
    metadata = merge(Dict(:sampler => :NestedSamplers, :sampler_type=>sampler),
                     Dict(Pair.(keys(stats), values(stats)))
                     )

    return PosteriorSamples(chain, (;weights); metadata)
end


end
