using .AdvancedHMC
export  HMC

Base.@kwdef struct HMC{S,I,P,T,A,D}
    metric::S
    integrator::I = Leapfrog
    trajectory::P = MultinomialTS
    termination::T = StrictGeneralisedNoUTurn(10, 1000.0)
    adaptor::A = StanHMCAdaptor
    targetacc::Float64 = 0.8
    autodiff::D = AD.ForwardDiffBackend()
end

samplertype(::Type{<:HMC}) = IsMCMC()

function AbstractMCMC.sample(tpost::TransformedPosterior, sampler::HMC, nsamples, args...;
                             init_params=nothing,
                             kwargs...)
    ℓ(x) = logdensity(tpost, x)

    function ∇ℓ(x)
        res = AD.value_and_gradient(sampler.autodiff, ℓ, x)
        return (first(res), first(last(res)))
    end

    θ0 = init_params
    if isnothing(init_params)
        @warn "No starting location chosen, picking start from random"
        θ0 = transform(tpost.transform, rand(tpost.prior))
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
    #return TupleVector(transform.(Ref(tpost), chain)), TupleVector(stats)
end
