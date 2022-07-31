module ComradeAHMC

using AbstractDifferentiation
using AbstractMCMC
using Reexport
@reexport using AdvancedHMC
using Comrade
using TupleVectors
using ArgCheck: @argcheck
using Random

export sample, AHMC

"""
    AdvancedHMC

Creates a sampler that uses the `AdvancedHMC` framework to construct
an Hamiltonian Monte Carlo NUTS sampler.

The user must specify the `metric` they want to use. Typically we recommend
`DiagEuclideanMetric` as a reasonable starting place. The other options
are chosen to match the `Stan` languages defaults and should provide
a good starting point. Please see the [AdvancedHMC docs](https://github.com/TuringLang/AdvancedHMC.jl)
for more information.

# Notes
For `autodiff` the user can either specify an `AbstractDifferentiation` backend. For
smaller geometric models we recommend `AD.ForwardDiffBackend` while for imaging we
recommend `AD.Zygote`.


# Fields

```julia
struct AHMC{S,I,P,T,A,D}
    `AdvancedHMC` metric to use
    metric::S
    `AdvancedHMC` integrator
    Defaults to `AdvancedHMC.Leapfrog`
    integrator::I = Leapfrog
    HMC trajectory sampler
    Defaults to `AdvancedHMC.MultinomialTS`
    trajectory::P = MultinomialTS
    HMC termination condition
    Defaults to `AdvancedHMC.StrictGeneralisedNoUTurn`
    termination::T = StrictGeneralisedNoUTurn(10, 1000.0)
    Adaptation strategy for mass matrix and stepsize
    Defaults to `AdvancedHMC.StanHMCAdaptor`
    adaptor::A = StanHMCAdaptor
    Target acceptance rate for all trajectories on the tree
    Defaults to 0.85
    targetacc::Float64 = 0.8
    `AbstractDifferentiation` autodiff backend
    Defaults to `AbstractDifferentiation.ForwardDiffBackend()`
    autodiff::D = AD.ForwardDiffBackend()
end

```
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
    `AbstractDifferentiation` autodiff backend
    Defaults to `AbstractDifferentiation.ForwardDiffBackend()`
    """
    autodiff::D = AD.ForwardDiffBackend()
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

This returns a tuple where the first element is `nchains` of `TupleVector`'s
each which contains the MCMC samples of one of the parallel chain and the second argument
is a set of ancilliary information about each set of samples.

# Notes

This will automatically transform the posterior to the flattened unconstrained space.
"""
function AbstractMCMC.sample(post::Comrade.Posterior,
    sampler::A, parallel::AbstractMCMC.AbstractMCMCEnsemble,
    nsamples, nchains;
    init_params=nothing, kwargs...
    ) where {A<:AHMC}
    tpost = asflat(post)
    if isnothing(init_params)
        θ0 = prior_sample(tpost, nchains)
    else
        θ0 = Comrade.HypercubeTransform.inverse.(Ref(tpost), init_params)
    end
    return Comrade.sample(tpost,
                            sampler,
                            parallel,
                            nsamples,
                            nchains;
                            init_params=θ0,
                            kwargs...)
end



function AbstractMCMC.sample(tpost::Comrade.TransformedPosterior,
                             sampler::AHMC, parallel::AbstractMCMC.AbstractMCMCEnsemble,
                             nsamples, nchains;
                             init_params=nothing, kwargs...
                             )

    ℓ = logdensityof(tpost)

    ∇ℓ = Comrade.make_pullback(ℓ, sampler.autodiff)
    θ0 = _initialize_hmc(tpost, init_params, nchains)
    model = AdvancedHMC.DifferentiableDensityModel(ℓ, ∇ℓ)
    metric = sampler.metric
    # This is a hack to get a good initial step size
    hamiltonian = Hamiltonian(metric, ℓ, ∇ℓ)
    ϵ0 = find_good_stepsize(hamiltonian, first(θ0))
    integrator = sampler.integrator(ϵ0)

    # form the HMCKernel
    kernel = HMCKernel(Trajectory{sampler.trajectory}(integrator, sampler.termination))
    adaptor = sampler.adaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(sampler.targetacc, integrator))

    res = AbstractMCMC.sample(Random.GLOBAL_RNG, model, kernel, metric, adaptor, parallel, nsamples, nchains; init_params=θ0, chain_type=Array, kwargs...)

    stats = [TupleVector(getproperty.(r, :stat)) for r in res]
    samples = [getproperty.(getproperty.(r, :z), :θ) for r in res]
    chains = [TupleVector(transform.(Ref(tpost), s)) for s in samples]
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

This returns a tuple where the first element is a `TupleVector` of the MCMC samples in parameter space
and the second argument is a set of ancilliary information about the sampler.


# Notes

This will automatically transform the posterior to the flattened unconstrained space.
"""
function AbstractMCMC.sample(tpost::Comrade.TransformedPosterior, sampler::AHMC, nsamples, args...;
                             init_params=nothing,
                             kwargs...)
    ℓ = logdensityof(tpost)

    ∇ℓ = Comrade.make_pullback(ℓ, sampler.autodiff)

    θ0 = init_params
    if isnothing(init_params)
        @warn "No starting location chosen, picking start from prior"
        θ0 = prior_sample(post)
    end
    model = AdvancedHMC.DifferentiableDensityModel(ℓ, ∇ℓ)
    metric = sampler.metric
    # This is a hack to get a good initial step size
    hamiltonian = Hamiltonian(metric, ℓ, ∇ℓ)
    ϵ0 = find_good_stepsize(hamiltonian, θ0)
    integrator = sampler.integrator(ϵ0)

    # form the HMCKernel
    kernel = HMCKernel(Trajectory{sampler.trajectory}(integrator, sampler.termination))
    adaptor = sampler.adaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(sampler.targetacc, integrator))

    res = AbstractMCMC.sample(model, kernel, metric, adaptor, nsamples, args...; init_params=θ0, chain_type=Array, kwargs...)

    stats = TupleVector(getproperty.(res, :stat))
    samples = getproperty.(getproperty.(res, :z), :θ)
    chain = transform.(Ref(tpost), samples)
    return TupleVector(chain), stats
end



end
