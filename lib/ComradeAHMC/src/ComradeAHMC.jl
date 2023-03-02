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

export sample, AHMC

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
    AbstractMCMC.sample(post::Comrade.Posterior,
                        sampler::AHMC,
                        nsamples;
                        init_params=nothing,
                        kwargs...)

Samples the posterior `post` using the AdvancedHMC sampler specified by `AHMC`.
This will run the sampler for `nsamples`.

To initialize the chain the user can set `init_params` to `Vector{NamedTuple}` whose
elements are the starting locations for each of the `nchains`. If no starting location
is specified `nchains` random samples from the prior will be chosen for the starting locations.

For possible `kwargs` please see the [`AdvancedHMC.jl docs`](https://github.com/TuringLang/AdvancedHMC.jl)

This returns a tuple where the first element is a `TypedTable` of the MCMC samples in parameter space
and the second argument is a set of ancilliary information about the sampler.


# Notes

This will automatically transform the posterior to the flattened unconstrained space.
"""
function AbstractMCMC.sample(rng::Random.AbstractRNG, tpost::Comrade.TransformedPosterior, sampler::AHMC, nsamples, args...;
                             init_params=nothing,
                             kwargs...)
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



end
