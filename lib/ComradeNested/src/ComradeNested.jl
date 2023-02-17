module ComradeNested

using Comrade
using AbstractMCMC
using InferenceObjects
using Reexport
using Random

@reexport using NestedSamplers

Comrade.samplertype(::Type{<:Nested}) = Comrade.IsCube()

"""
    AbstractMCMC.sample(post::Comrade.Posterior, smplr::Nested, args...; kwargs...)

Sample the posterior `post` using `NestedSamplers.jl` `Nested` sampler. The `args/kwargs`
are forwarded to `NestedSampler` for more information see its [docs](https://github.com/TuringLang/NestedSamplers.jl)

This returns an InferenceData object which includes the weighted samples in the posterior group
and the weights of those samples in the sample_stats group.

To create equally weighted samples the user can use
```julia
samples = sample(post, Nested(dimension(post), 1000))
equal_weighted_chain = ComradeNested.equalresample(samples, 1000)
```

"""
function AbstractMCMC.sample(rng::Random.AbstractRNG, post::Comrade.TransformedPosterior, sampler::Nested, args...; kwargs...)
    ℓ = logdensityof(post)
    model = NestedModel(ℓ, identity)

    samples, stats = sample(rng, model, sampler, args...; chain_type=Array, kwargs...)

    sample_stats = (weights = samples[:, end],)
    attrs = Dict(
        :logz    => stats.logz,
        :logzerr => stats.logzerr,
        :logl    => stats.logl,
        )
    chain = transform.(Ref(post), eachrow(@view samples[:,1:end-1]))
    library = "Comrade NestedSamplers.jl"
    return from_namedtuple([chain]; sample_stats, library, attrs)
end

function equalresample(res::InferenceData, nsamples)
    weights = vec(res.sample_stats.weights)
    eq = sample(Tables.rowtable(res.posterior[chain=1]), Weights(weights), nsamples)
    return from_namedtuple([eq])
end


end
