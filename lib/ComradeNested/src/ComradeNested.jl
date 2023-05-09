module ComradeNested

using Pyehtim
using Comrade
using AbstractMCMC
using TypedTables
using Reexport
using Random

@reexport using NestedSamplers

Comrade.samplertype(::Type{<:Nested}) = Comrade.IsCube()

"""
    AbstractMCMC.sample(post::Comrade.Posterior, smplr::Nested, args...; kwargs...)

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
function AbstractMCMC.sample(rng::Random.AbstractRNG, post::Comrade.TransformedPosterior, sampler::Nested, args...; kwargs...)
    ℓ = logdensityof(post)
    model = NestedModel(ℓ, identity)

    samples, stats = sample(rng, model, sampler, args...; chain_type=Array, kwargs...)
    weights = samples[:, end]
    chain = transform.(Ref(post), eachrow(samples[:,1:end-1]))
    return Table(chain), merge((;weights,), stats)
end

end
