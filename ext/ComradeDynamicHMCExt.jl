module ComradeDynamicHMCExt

using Comrade
using DynamicHMC

using DocStringExtensions
using HypercubeTransform
using Random
using LogDensityProblems
using Serialization
using StatsBase


function DynamicHMC.mcmc_with_warmup(
        rng, post::Comrade.VLBIPosterior, N;
        kwargs...
    )

    if isnothing(Comrade.admode(post))
        throw(ArgumentError("You must specify an automatic differentiation type in VLBIPosterior with admode kwarg"))
    else
        tpost = asflat(post)
    end

    results = mcmc_with_warmup(rng, tpost, N; kwargs...)

    return PosteriorSamples(
        transform.(Ref(tpost), eachcol(results.posterior_matrix)),
        results.tree_statistics;
        metadata = Dict(
            :sampler = :DynamicHMC, :ϵ => results.ϵ,
            :mass_matrix => results.κ
        )
    )
end

function DynamicHMC.mcmc_with_warmup(
        rng, post::Comrade.VLBIPosterior, N, output::DiskStore;
        kwargs...
    )

    (; name, stride) = output
    stride = min(stride, N)
    nscans = nsamples ÷ output_stride + (nsamples % output_stride != 0 ? 1 : 0)
    outbase = joinpath(name, "samples", "output_scan_")


    tpost = asflat(post)
    results = mcmc_with_warmup(rng, tpost, N; kwargs...)

    return PosteriorSamples(
        transform.(Ref(tpost), eachcol(results.posterior_matrix)),
        results.tree_statistics;
        metadata = Dict(
            :sampler = :DynamicHMC, :ϵ => results.ϵ,
            :mass_matrix => results.κ
        )
    )
end


end
