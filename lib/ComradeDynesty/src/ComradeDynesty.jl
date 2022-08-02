module ComradeDynesty

using Comrade

using AbstractMCMC
using TypedTables
using Reexport

@reexport using Dynesty


Comrade.samplertype(::Type{<:NestedSampler}) = Comrade.IsCube()
Comrade.samplertype(::Type{<:DynamicNestedSampler}) = Comrade.IsCube()

"""
    AbstractMCMC.sample(post::Comrade.Posterior, smplr::Dynesty.NestedSampler, args...; kwargs...)
    AbstractMCMC.sample(post::Comrade.Posterior, smplr::Dynesty.DynamicNestedSampler, args...; kwargs...)

Sample the posterior `post` using `Dynesty.jl` `NestedSampler/DynamicNestedSampler` sampler.
The `args/kwargs`
are forwarded to `Dynesty` for more information see its [docs](https://github.com/ptiede/Dynesty.jl)

This returns a tuple where the first element are the weighted samples from dynesty in a TypedTable.
The second element includes additional information about the samples, like the log-likelihood,
evidence, evidence error, and the sample weights.

To create equally weighted samples the user can use
```julia
using StatsBase
chain, stats = sample(post, NestedSampler(dimension(post), 1000))
equal_weighted_chain = sample(chain, Weights(stats.weights), 10_000)
```
"""
function AbstractMCMC.sample(post::Comrade.TransformedPosterior,
                             sampler::Union{NestedSampler, DynamicNestedSampler},
                             args...;
                             kwargs...)
    ℓ = logdensityof(post)
    kw = delete!(Dict(kwargs), :init_params)
    res = sample(ℓ, identity, sampler, args...; kw...)
    samples, weights = res["samples"].T, exp.(res["logwt"].T .- res["logz"][end])
    chain = transform.(Ref(post), eachcol(samples)) |> Table
    stats = (logl = res["logl"].T,
             logz = res["logz"][end],
             logzerr = res["logz"][end],
             weights = weights,
            )
    return Table(chain), stats
end


end
