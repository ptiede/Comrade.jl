module ComradeDynesty

using Comrade
using Dynesty

using AbstractMCMC
using Reexport
using Random

@reexport using Dynesty


Comrade.samplertype(::Type{<:NestedSampler}) = Comrade.IsCube()
Comrade.samplertype(::Type{<:DynamicNestedSampler}) = Comrade.IsCube()

"""
    AbstractMCMC.sample(post::Comrade.VLBIPosterior, smplr::Dynesty.NestedSampler, args...; kwargs...)
    AbstractMCMC.sample(post::Comrade.VLBIPosterior, smplr::Dynesty.DynamicNestedSampler, args...; kwargs...)

Sample the posterior `post` using `Dynesty.jl` `NestedSampler/DynamicNestedSampler` sampler.
The `args/kwargs`
are forwarded to `Dynesty` for more information see its [docs](https://github.com/ptiede/Dynesty.jl)

This returns a PosteriorSamples object.
The `samplerstats` includes additional information about the samples, like the log-likelihood,
evidence, evidence error, and the sample weights. The final element of the tuple is the original
dynesty output file.

To create equally weighted samples the user can use
```julia
using StatsBase
chain = sample(post, NestedSampler(dimension(post), 1000))
equal_weighted_chain = sample(chain, Weights(stats.weights), 10_000)
```
"""
function AbstractMCMC.sample(::Random.AbstractRNG, post::Comrade.TransformedVLBIPosterior,
                             sampler::Union{NestedSampler, DynamicNestedSampler}
                             ; initial_params=nothing,
                             kwargs...)
    ℓ = logdensityof(post)
    kw = delete!(Dict(kwargs), :initial_params)
    res = dysample(ℓ, identity, sampler; kw...)
    # Make sure that res["sample"] is an array and use transpose
    samples, weights = transpose(Dynesty.PythonCall.pyconvert(Array, res["samples"])),
                       exp.(Dynesty.PythonCall.pyconvert(Vector, res["logwt"] - res["logz"][-1]))
    chain = transform.(Ref(post), eachcol(samples))
    stats = (logl = Dynesty.PythonCall.pyconvert(Vector, res["logl"]),
             weights = weights,
            )

    logz = Dynesty.PythonCall.pyconvert(Float64, res["logz"][-1])
    logzerr = Dynesty.PythonCall.pyconvert(Float64, res["logzerr"][-1])

    return PosteriorSamples(chain, stats; metadata=Dict(:sampler => :dynesty, :dynesty_output => res, :logz => logz, :logzerr => logzerr))
end



end
