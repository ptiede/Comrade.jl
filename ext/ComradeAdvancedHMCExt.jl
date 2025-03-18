module ComradeAdvancedHMCExt

using Comrade
using AdvancedHMC
using AdvancedHMC: AbstractHMCSampler

using AbstractMCMC
using AbstractMCMC: Sample
using Accessors
using ArgCheck
using DocStringExtensions
using HypercubeTransform
using LogDensityProblems
using Printf
using Random
using StatsBase
using Serialization


function initialize_params(tpost, initial_params)
    isnothing(initial_params) && return prior_sample(tpost)
    return HypercubeTransform.inverse(tpost, initial_params)
end

# internal method that makes the HMCSampler
function make_sampler(rng, ∇ℓ, sampler::AdvancedHMC.AbstractHMCSampler, θ0)
    metric = AdvancedHMC.make_metric(sampler, ∇ℓ)
    hamil  = AdvancedHMC.Hamiltonian(metric, ∇ℓ)
    ϵ      = AdvancedHMC.make_step_size(rng, sampler, hamil, θ0)
    integr = AdvancedHMC.make_integrator(sampler, ϵ)
    κ      = AdvancedHMC.make_kernel(sampler, integr)
    adaptor= AdvancedHMC.make_adaptor(sampler, metric, integr)
    return AdvancedHMC.LogDensityModel(∇ℓ), HMCSampler(κ, metric, adaptor)
end


function AbstractMCMC.Sample(
            rng::Random.AbstractRNG, tpost::Comrade.TransformedVLBIPosterior,
            sampler::AbstractHMCSampler; initial_params=nothiing, kwargs...)
    θ0 = initialize_params(tpost, initial_params)
    model, smplr = make_sampler(rng, tpost, sampler, θ0)
    return AbstractMCMC.Sample(rng, model, smplr; initial_params=θ0, kwargs...)
end

"""
    sample(rng, post::VLBIPosterior, sampler::AbstractHMCSampler, nsamples, args...;saveto=MemoryStore(), initial_params=nothing, kwargs...)

Sample from the posterior `post` using the sampler `sampler` for `nsamples` samples. Additional
arguments are forwarded to AbstractMCMC.sample. If `saveto` is a DiskStore, the samples will be
saved to disk. If `initial_params` is not `nothing` then the sampler will start from that point.

## Arguments

 - rng: The random number generator to use
 - post: The posterior to sample from
 - nsamples: The total number of samples

## Keyword Arguments

 - `saveto`: If a DiskStore, the samples will be saved to disk, if [`MemoryStore`](@ref) the samples will be stored in memory/ram.
 - `initial_params`: The initial parameters to start the sampler from. If `nothing` then the sampler will start from a random point in the prior.
 - `kwargs`: Additional keyword arguments to pass to the sampler. Examples include `n_adapts` which is the total number of samples to use for adaptation.
    To see the others see the AdvancedHMC documentation.
"""
function AbstractMCMC.sample(
        rng::Random.AbstractRNG, post::Comrade.VLBIPosterior,
        sampler::AbstractHMCSampler, nsamples, args...;
        saveto=MemoryStore(), initial_params=nothing, kwargs...)

    saveto isa DiskStore && return sample_to_disk(rng, post, sampler, nsamples, args...; outdir=saveto.name, output_stride=min(saveto.stride, nsamples), initial_params, kwargs...)

    if isnothing(Comrade.admode(post))
        throw(ArgumentError("You must specify an automatic differentiation type in VLBIPosterior with admode kwarg"))
    else
        tpost = asflat(post)
    end

    tpost = asflat(post)
    θ0 = initialize_params(tpost, initial_params)
    model, smplr = make_sampler(rng, tpost, sampler, θ0)

    res = sample(rng, model, smplr, nsamples, args...;
                 initial_params=θ0, saveto=saveto, chain_type=Array, kwargs...)

    stats = getproperty.(res, :stat)
    samples = getproperty.(getproperty.(res, :z), :θ)
    chain = transform.(Ref(tpost), samples)
    return PosteriorSamples(chain, stats; metadata=Dict(:sampler=>:AHMC, :sampler_kwargs=>kwargs, :post=>tpost))
end

# Disk sampling stuff goes here
function initialize(rng::Random.AbstractRNG, tpost::Comrade.TransformedVLBIPosterior,
    sampler::AbstractHMCSampler, nsamples, outbase, args...;
    n_adapts = min(nsamples÷2, 1000),
    initial_params=nothing, outdir = "Results",
    output_stride=min(100, nsamples),
    restart = false,
    kwargs...)

    if restart
        @assert isfile(joinpath(outdir, "checkpoint.jls")) "Checkpoint file does not exist in $(outdir)"
        tmp = deserialize(joinpath(outdir, "checkpoint.jls"))
        (;pt, state, out, iter) = tmp
        if iter*output_stride >= nsamples
            @warn("Not sampling because the current number of samples is greater than the number you requested")
            return pt, state, out, iter
        end
        if pt.c.coll.stop != nsamples
            @warn("The number of samples wanted in the stored checkpoint does not match the number of samples requested."*
                  "Changing to the number you requested")
            pt = @set pt.c.coll.stop = nsamples
        end
        @info "Resuming from checkpoint on iteration $iter"
        return pt, state, out, iter
    end

    mkpath(joinpath(outdir, "samples"))
    θ0 = initial_params
    if isnothing(initial_params)
        @warn "No starting location chosen, picking start from prior"
        θ0 = prior_sample(rng, tpost)
    end
    t = Sample(rng, tpost, sampler; initial_params=θ0, n_adapts, kwargs...)(1:nsamples)
    pt = Iterators.partition(t, output_stride)
    nscans = nsamples÷output_stride + (nsamples%output_stride!=0 ? 1 : 0)

    # Now save the output
    out = Comrade.DiskOutput(abspath(outdir), nscans, output_stride, nsamples)
    serialize(joinpath(outdir, "parameters.jls"), (;params=out))

    tmp = @timed iterate(pt)
    state, iter = _process_samples(pt, tpost, tmp.value, tmp.time, nscans, out, outbase, outdir, 1)
    return pt, state, out, iter
end



function _process_samples(pt, tpost, next, time, nscans, out, outbase, outdir, iter)
    (chain, state) = next
    stats = getproperty.(chain, :stat)
    samples = transform.(Ref(tpost), getproperty.(getproperty.(chain, :z), :θ))
    s = PosteriorSamples(samples, stats; metadata=Dict(:sampler=>:AHMC))
    @info("Scan $(iter)/$nscans statistics:\n"*
          "  Total time:     $(time) seconds\n"*
          "  Mean tree depth: $(round(mean(samplerstats(s).tree_depth); digits=1))\n"*
          "  Mode tree depth: $(round(StatsBase.mode(samplerstats(s).tree_depth); digits=1))\n"*
          "  n-divergences:   $(sum(samplerstats(s).numerical_error))/$(length(stats))\n"*
          "  Avg log-post:    $(mean(samplerstats(s).log_density))\n")

    serialize(outbase*(@sprintf "%05d.jls" iter), (samples=Comrade.postsamples(s), stats=Comrade.samplerstats(s)))
    # chain = nothing
    # samples = nothing
    # stats = nothing
    # s = nothing
    # GC.gc()
    iter += 1
    serialize(joinpath(outdir, "checkpoint.jls"), (;pt, state, out, iter))
    return state, iter
end

function sample_to_disk(rng::Random.AbstractRNG, post::Comrade.VLBIPosterior,
                        sampler::AbstractHMCSampler, nsamples, args...;
                        n_adapts = min(nsamples÷2, 1000),
                        initial_params=nothing, outdir = "Results",
                        restart=false,
                        output_stride=min(100, nsamples),
                        kwargs...)


    tpost = asflat(post)
    nscans = nsamples÷output_stride + (nsamples%output_stride!=0 ? 1 : 0)
    outbase = joinpath(outdir, "samples", "output_scan_")

    pt, state, out, i = initialize(
                            rng, tpost, sampler, nsamples, outbase, args...;
                            n_adapts,
                            initial_params, restart, outdir, output_stride, kwargs...
                        )

    tmp = @timed iterate(pt, state)
    t = tmp.time
    next = tmp.value
    while !isnothing(next)
        t = @elapsed begin
            state, i = _process_samples(pt, tpost, next, t, nscans, out, outbase, outdir, i)
            next = iterate(pt, state)
        end
    end


    return out
end





end
