module ComradeAHMC

using AbstractMCMC
using Reexport
@reexport using AdvancedHMC
using Comrade
using DocStringExtensions
using LogDensityProblems, LogDensityProblemsAD
using TypedTables
using ArgCheck: @argcheck
using Random
using JLD2
using Printf
using AbstractMCMC: Sample
export sample, AHMC, Sample, Memory, Disk, load_table

"""
    AHMC

Creates a sampler that uses the `AdvancedHMC` framework to construct
an Hamiltonian Monte Carlo NUTS sampler.

The user must specify the `metric` they want to use. Typically we recommend
`DiagEuclideanMetric` as a reasonable starting place. The other options
are chosen to match the `Stan` languages defaults and should provide
a good starting point. Please see the [AdvancedHMC docs](https://github.com/TuringLang/AdvancedHMC.jl)
for more information.

# Notes
For `autodiff` the must provide a `Val(::Symbol)` that specifies the AD backend. Currently,
we use `LogDensityProblemsAD`.


# Fields
$(FIELDS)

"""
Base.@kwdef struct AHMC{S,I,P,T,A,D}
    """
    `AdvancedHMC` metric to use
    """
    metric::S
    """
    `AdvancedHMC` integrator
    Defaults to `AdvancedHMC.Leapfrog`
    """
    integrator::I = Leapfrog
    """
    HMC trajectory sampler
    Defaults to `AdvancedHMC.MultinomialTS`
    """
    trajectory::P = MultinomialTS
    """
    HMC termination condition
    Defaults to `AdvancedHMC.StrictGeneralisedNoUTurn`
    """
    termination::T = StrictGeneralisedNoUTurn(10, 1000.0)
    """
    Adaptation strategy for mass matrix and stepsize
    Defaults to `AdvancedHMC.StanHMCAdaptor`
    """
    adaptor::A = StanHMCAdaptor
    """
    Target acceptance rate for all trajectories on the tree
    Defaults to 0.85
    """
    targetacc::Float64 = 0.8
    """
    The number of steps for the initial tuning phase.
    Defaults to 75 which is the Stan default
    """
    init_buffer::Int = 75
    """
    The number of steps for the final fast step size adaptation
    Default if 50 which is the Stan default
    """
    term_buffer::Int = 50
    """
    The number of steps to tune the covariance before the first doubling
    Default is 23 which is the Stan default
    """
    window_size::Int = 25
    """
    autodiff backend see [`LogDensitProblemsAD.jl`](https://github.com/tpapp/LogDensityProblemsAD.jl)
    for possible backends. The default is `Zygote` which is appropriate for high dimensional problems.
    """
    autodiff::D = Val(:Zygote)
end


Comrade.samplertype(::Type{<:AHMC}) = Comrade.IsFlat()

function _initialize_hmc(tpost::Comrade.TransformedPosterior, init_params, nchains)
    isnothing(init_params) && return Comrade.HypercubeTransform.inverse.(Ref(tpost.transform), rand(tpost.lpost.prior, nchains))
    @argcheck length(init_params) == nchains
    return init_params
end

"""
    AbstractMCMC.sample(post::Comrade.Posterior,
                        sampler::AHMC,
                        parallel::AbstractMCMC.AbstractMCMCEnsemble,
                        nsamples,
                        nchainsl;
                        init_params=nothing,
                        kwargs...)

Samples the posterior `post` using the AdvancedHMC sampler specified by `AHMC`.
This will sample `nchains` copies of the posterior using the `parallel` scheme.
Each chain will be sampled for `nsamples`.

To initialize the chain the user can set `init_params` to `Vector{NamedTuple}` whose
elements are the starting locations for each of the `nchains`. If no starting location
is specified `nchains` random samples from the prior will be chosen for the starting locations.

For possible `kwargs` please see the [`AdvancedHMC.jl docs`](https://github.com/TuringLang/AdvancedHMC.jl)

This returns a tuple where the first element is `nchains` of `TypedTable`'s
each which contains the MCMC samples of one of the parallel chain and the second argument
is a set of ancilliary information about each set of samples.

# Notes

This will automatically transform the posterior to the flattened unconstrained space.
"""
function AbstractMCMC.sample(
    rng::Random.AbstractRNG, post::Comrade.Posterior,
    sampler::A, parallel::AbstractMCMC.AbstractMCMCEnsemble,
    nsamples, nchains;
    init_params=nothing, kwargs...
    ) where {A<:AHMC}
    tpost = asflat(post)
    if isnothing(init_params)
        θ0 = prior_sample(rng, tpost, nchains)
    else
        θ0 = Comrade.HypercubeTransform.inverse.(Ref(tpost), init_params)
    end
    return Comrade.sample(rng,
                            tpost,
                            sampler,
                            parallel,
                            nsamples,
                            nchains;
                            init_params=θ0,
                            kwargs...)
end

function make_sampler(∇ℓ, sampler::AHMC, θ0)
    model = AdvancedHMC.LogDensityModel(∇ℓ)
    initial_ϵ = 1e-4
    integrator = sampler.integrator(initial_ϵ)
    proposal = HMCKernel(Trajectory{sampler.trajectory}(integrator, sampler.termination))
    adaptor = sampler.adaptor(MassMatrixAdaptor(sampler.metric), StepSizeAdaptor(sampler.targetacc, integrator);
                        init_buffer = sampler.init_buffer,
                        term_buffer = sampler.term_buffer,
                        window_size = sampler.window_size)

    return model, proposal, sampler.metric, adaptor
end

function AbstractMCMC.Sample(rng::Random.AbstractRNG, tpost::Comrade.TransformedPosterior, sampler::AHMC; kwargs...)
    ∇ℓ = ADgradient(sampler.autodiff, tpost)
    model, proposal, metric, adaptor = make_sampler(∇ℓ, sampler, 0)
    return AbstractMCMC.Sample(rng, model, AdvancedHMC.HMCSampler(proposal, metric, adaptor); kwargs...)
end


function AbstractMCMC.sample(rng::Random.AbstractRNG, tpost::Comrade.TransformedPosterior,
                             sampler::AHMC, parallel::AbstractMCMC.AbstractMCMCEnsemble,
                             nsamples, nchains;
                             init_params=nothing, kwargs...
                             )


    ∇ℓ = ADgradient(sampler.autodiff, tpost)
    θ0 = _initialize_hmc(tpost, init_params, nchains)
    model, proposal, metric, adaptor = make_sampler(∇ℓ, sampler, first(θ0))


    res = AbstractMCMC.sample(
                rng,
                model, proposal,
                metric, adaptor,
                parallel, nsamples, nchains;
                init_params=θ0,
                chain_type = Array, kwargs...
                )

    stats = [Table(getproperty.(r, :stat)) for r in res]
    samples = [getproperty.(getproperty.(r, :z), :θ) for r in res]
    chains = [Table(transform.(Ref(tpost), s)) for s in samples]
    return chains, stats

end

"""
    Memory

Stored the HMC samplers in memory or ram.
"""
struct Memory end


"""
    Disk

Type that specifies to save the HMC results to disk.

# Fields
$(FIELDS)
"""
Base.@kwdef struct Disk
    """
    Path of the directory where the results will be saved. If the path does not exist
    it will be automatically created.
    """
    name::String = "Results"
    """
    The output stride, i.e. every `stride` steps the MCMC output will be dumped to disk.
    """
    stride::Int = 500
end
Disk(name::String) = Disk(name, 500)

"""
    AbstractMCMC.sample(post::Comrade.Posterior,
                        sampler::AHMC,
                        nsamples;
                        init_params=nothing,
                        saveto::Union{Memory, Disk}=Memory(),
                        kwargs...)

Samples the posterior `post` using the AdvancedHMC sampler specified by `AHMC`.
This will run the sampler for `nsamples`.

To initialize the chain the user can set `init_params` to `Vector{NamedTuple}` whose
elements are the starting locations for each of the `nchains`. If no starting location
is specified `nchains` random samples from the prior will be chosen for the starting locations.

With `saveto` the user can optionally specify whether to store the samples in memory `Memory`
or save directly to disk with `Disk(filename, stride)`. The `stride` controls how often t
he samples are dumped to disk.

For possible `kwargs` please see the [`AdvancedHMC.jl docs`](https://github.com/TuringLang/AdvancedHMC.jl)

This returns a tuple where the first element is a `TypedTable` of the MCMC samples in parameter space
and the second argument is a set of ancilliary information about the sampler.


# Notes

This will automatically transform the posterior to the flattened unconstrained space.
"""
function AbstractMCMC.sample(rng::Random.AbstractRNG, tpost::Comrade.TransformedPosterior,
                             sampler::AHMC, nsamples, args...;
                             saveto=Memory(),
                             init_params=nothing,
                             kwargs...)


    saveto isa Disk && return sample_to_disk(rng, tpost, sampler, nsamples, args...; outdir=saveto.name, output_stride=min(saveto.stride, nsamples), init_params, kwargs...)

    ℓ = logdensityof(tpost)

    ∇ℓ = ADgradient(sampler.autodiff, tpost)
    θ0 = init_params

    if isnothing(init_params)
        @warn "No starting location chosen, picking start from prior"
        θ0 = prior_sample(rng, post)
    end

    model, proposal, metric, adaptor = make_sampler(∇ℓ, sampler, first(θ0))


    res = AbstractMCMC.sample(
                rng,
                model, proposal,
                metric, adaptor,
                nsamples;
                init_params=θ0,
                chain_type = Array, kwargs...
                )

    stats = Table(getproperty.(res, :stat))
    samples = getproperty.(getproperty.(res, :z), :θ)
    chain = transform.(Ref(tpost), samples)
    return Table(chain), stats
end

struct DiskOutput
    filename::String
    nfiles::Int
    stride::Int
    nsamples::Int
end


function sample_to_disk(rng::Random.AbstractRNG, tpost::Comrade.TransformedPosterior,
                        sampler::AHMC, nsamples, args...;
                        init_params=nothing, outdir = "Results",
                        output_stride=min(100, nsamples),
                        kwargs...)


    mkpath(outdir)
    θ0 = init_params
    if isnothing(init_params)
        @warn "No starting location chosen, picking start from prior"
        θ0 = prior_sample(rng, tpost)
    end
    t = Sample(rng, tpost, sampler; init_params, kwargs...)(1:nsamples)
    pt = Iterators.partition(t, output_stride)
    outbase = joinpath(outdir, "samples", "output_scan_")
    nscans = nsamples÷output_stride + (nsamples%output_stride!=0 ? 1 : 0)

    next = iterate(pt)
    i = 1
    while !isnothing(next)
        (chain, state) = next
        t = @elapsed begin
            stats = Table(getproperty.(chain, :stat))
            samples = transform.(Ref(tpost), getproperty.(getproperty.(chain, :z), :θ)) |> Table
            jldsave(outbase*(@sprintf "%05d.jld2" i); stats, samples)
            next = iterate(pt, state)
        end
        @info "On scan $i/$nscans it took $(t) seconds"
        i += 1
    end

    # Now save the output as well
    out = DiskOutput(outdir, nscans, output_stride, nsamples)
    jldsave(joinpath(outdir, "parameters"), out)

    return out
end

"""
    load_table(out::DiskOutput, indices::Union{Base.Colon, UnitRange, StepRange}=Base.Colon(); table="samples")
    load_table(out::String, indices::Union{Base.Colon, UnitRange, StepRange}=Base.Colon(); table="samples")

The the results from a HMC run saved to disk. To read in the output the user can either
pass the resulting `out` object, or the path to the directory that the results were saved,
i.e. the path specified in [`Disk`](@ref).

By default if just a `out` object is given the entire table of samples will be read in.
Otherwise the use can specify the specific range
of `indices` they are interested in. Finally, HMC save two tables, `samples` and `stats`. The
`samples` table correspond to the actual MCMC samples from HMC stored as a `TypedTable`.
The `stats` table hold the rest of the sampling statistics such as the evaluated log-density,
the number of leapfrog steps, etc.
"""
function load_table(
        out::DiskOutput,
        indices::Union{Base.Colon, UnitRange, StepRange}=Base.Colon(); table="samples"
        )
    @assert (table == "samples" || table == "stats") "Please select either `samples` or `stats`"
    d = readdir(joinpath(out.filename, "samples"), join=true)

    # load and return the entire table
    indices == Base.Colon() && (return reduce(vcat, load.(d, table)))

    # Now get the index of the first file
    ind0 = first(indices)
    (ind0 < 1) && throw(BoundsError(1:out.nsamples, ind0))
    # Now let's find the file
    find0 = ind0÷out.stride + 1
    offset0 = ind0%out.stride # now get the offset
    if offset0 == 0
        find0 = find0 - 1
        offset0 = out.stride
    end

    ind1 = last(indices)
    (ind1 > out.nsamples) && throw(BoundsError(1:out.nsamples, ind1))
    find1 = ind1÷out.stride + 1
    offset1 = ind1%out.stride # now get the offset
    if offset1 == 0
        find1 = find1 - 1
        offset1 = out.stride
    end

    t = reduce(vcat, load.(d[find0:find1], table))
    return t[offset0:step(indices):(out.stride*(find1-find0) + offset1)]
end


end
