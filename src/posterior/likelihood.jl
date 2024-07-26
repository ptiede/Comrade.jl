struct ConditionedLikelihood{F, O}
    kernel::F
    obs::O
end
@inline DensityInterface.logdensityof(d::ConditionedLikelihood, μ) = logdensityof(@inline(d.kernel(μ)), d.obs)


"""
    likelihood(d::ConditionedLikelihood, μ)

Returns the likelihood of the model, with parameters μ. That is, we return the
distribution of the data given the model parameters μ. Samples from this distribution
are simulated data.
"""
likelihood(d::ConditionedLikelihood, μ) = d.kernel(μ)



# internal function that creates the likelihood for a set of complex visibilities
function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTVisibilityDatum})
    Σ = noise(data).^2
    vis = measurement(data)
    lnorm = VLBILikelihoods.lognorm(ComplexVisLikelihood(vis, Σ))
    ℓ = ConditionedLikelihood(vis) do μ
        ComplexVisLikelihood(μ, Σ, lnorm)
    end
    return ℓ
end

function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTCoherencyDatum})
    Σ = map(x->x.^2, noise(data))
    vis = measurement(data)
    # lnorm = VLBILikelihoods.lognorm(CoherencyLikelihood(vis, Σ))
    ℓ = ConditionedLikelihood(vis) do μ
        # @info typeof(μ)
        CoherencyLikelihood(baseimage(μ), Σ, 0.0)
    end
    return ℓ
end

# internal function that creates the likelihood for a set of visibility amplitudes
function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTVisibilityAmplitudeDatum})
    Σ = noise(data).^2
    amp = measurement(data)
    ℓ = ConditionedLikelihood(amp) do μ
        RiceAmplitudeLikelihood(abs.(μ), Σ)
    end
    return ℓ
end

# internal function that creates the likelihood for a set of log closure amplitudes
function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTLogClosureAmplitudeDatum})
    Σlca = factornoisecovariance(arrayconfig(data))
    f = Base.Fix2(logclosure_amplitudes, designmat(arrayconfig(data)))
    amp = measurement(data)
    lnorm = VLBILikelihoods.lognorm(AmplitudeLikelihood(amp, Σlca))
    ℓ = ConditionedLikelihood(amp) do μ
        AmplitudeLikelihood(f(μ), Σlca, lnorm)
    end
    return ℓ
end


# internal function that creates the likelihood for a set of closure phase datum
function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTClosurePhaseDatum})
    Σcp = factornoisecovariance(arrayconfig(data))
    f = Base.Fix2(closure_phases, designmat(arrayconfig(data)))
    phase = measurement(data)
    lnorm = VLBILikelihoods.lognorm(ClosurePhaseLikelihood(phase, Σcp))
    ℓ = ConditionedLikelihood(phase) do μ
        ClosurePhaseLikelihood(f(μ), Σcp, lnorm)
    end

    return ℓ
end
