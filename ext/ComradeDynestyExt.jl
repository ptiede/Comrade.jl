module ComradeDynestyExt

using AbstractMCMC
using Comrade
using Dynesty
using Random

"""
    dysample(post::Comrade.VLBIPosterior, smplr::Dynesty.NestedSampler; kwargs...)
    dysample(post::Comrade.VLBIPosterior, smplr::Dynesty.DynamicNestedSampler; kwargs...)

Sample the posterior `post` using `Dynesty.jl` `NestedSampler/DynamicNestedSampler` sampler.
The `kwargs` are forwarded to `Dynesty`; for more information see its
[docs](https://github.com/ptiede/Dynesty.jl).

This returns a `PosteriorSamples` object. The `samplerstats` hold the log-likelihood and the
weight of each sample. The `samplerinfo` holds the evidence `:logz`, its error `:logzerr`, and
the original dynesty output `:dynesty_output`.

To create equally weighted samples use [`resample_equal`](@ref):
```julia
chain = dysample(post, NestedSampler(; nlive = 1000))
equal_weighted_chain = resample_equal(chain, 10_000)
```
"""
function Dynesty.dysample(
        post::Comrade.VLBIPosterior,
        sampler::Union{Dynesty.NestedSampler, DynamicNestedSampler};
        kwargs...
    )

    tpost = ascube(post)
    ℓ = logdensityof(tpost)
    function l(x)
        out = ℓ(x)
        return isfinite(out) ? out : -1.0e300
    end
    res = dysample(l, identity, Comrade.dimension(tpost), sampler; kwargs...)
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
