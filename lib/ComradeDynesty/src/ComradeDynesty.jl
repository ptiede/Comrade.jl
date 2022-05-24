module ComradeDynesty

using Comrade

using AbstractMCMC
using TupleVectors
using Reexport

@reexport using Dynesty


Comrade.samplertype(::Type{<:NestedSampler}) = Comrade.IsCube()

function AbstractMCMC.sample(post::Comrade.TransformedPosterior, sampler::NestedSampler, args...; kwargs...)
    ℓ(x) = logdensityof(post, x)
    kw = delete!(Dict(kwargs), :init_params)
    res = sample(ℓ, identity, sampler, args...; kw...)
    samples, weights = res["samples"].T, exp.(res["logwt"].T .- res["logz"][end])
    chain = transform.(Ref(post), eachcol(samples)) |> TupleVector
    stats = (logl = res["logl"].T,
             logz = res["logz"][end],
             logzerr = res["logz"][end],
             weights = weights,
            )
    return chain, stats
end


end
