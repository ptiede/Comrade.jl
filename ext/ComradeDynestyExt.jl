module ComradeDynestyExt

using AbstractMCMC
using Comrade
using Dynesty
using Random

"""
    dysample(post::Comrade.VLBIPosterior, smplr::Dynesty.NestedSampler, args...; kwargs...)
    dysample(post::Comrade.VLBIPosterior, smplr::Dynesty.DynamicNestedSampler, args...; kwargs...)

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
equal_weighted_chain = sample(chain, Weights(samplerstats(chain).weights), 10_000)
```
"""
function Dynesty.dysample(
        post::Comrade.VLBIPosterior,
        sampler::Union{Dynesty.NestedSampler, DynamicNestedSampler};
        kwargs...
    )

    tpost = ascube(post)
    ℓ = logdensityof(tpost)
    res = dysample(ℓ, identity, Comrade.dimension(tpost), sampler; kwargs...)
    # Make sure that res["sample"] is an array and use transpose
    samples, weights = transpose(Dynesty.PythonCall.pyconvert(Array, res["samples"])),
        exp.(Dynesty.PythonCall.pyconvert(Vector, res["logwt"] - res["logz"][-1]))
    chain = transform.(Ref(tpost), eachcol(samples))
    stats = (
        logl = Dynesty.PythonCall.pyconvert(Vector, res["logl"]),
        weights = weights,
    )

    logz = Dynesty.PythonCall.pyconvert(Float64, res["logz"][-1])
    logzerr = Dynesty.PythonCall.pyconvert(Float64, res["logzerr"][-1])

    return PosteriorSamples(chain, stats; metadata = Dict(:sampler => :dynesty, :dynesty_output => res, :logz => logz, :logzerr => logzerr))
end


end
