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


struct _Visibility{S, L}
    S::S
    L::L
end

function (v::_Visibility)(μ)
    return ComplexVisLikelihood(baseimage(μ), v.S, v.L)
end

# internal function that creates the likelihood for a set of complex visibilities
function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTVisibilityDatum})
    Σ = noise(data) .^ 2
    vis = measurement(data)
    lnorm = VLBILikelihoods.lognorm(ComplexVisLikelihood(vis, Σ))
    ℓ = ConditionedLikelihood(_Visibility(Σ, lnorm), vis)
    return ℓ
end

struct _Coherency{S, L}
    S::S
    L::L
end

function (c::_Coherency)(μ)
    return CoherencyLikelihood(baseimage(μ), c.S, c.L)
end

function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTCoherencyDatum})
    Σ = map(x -> x .^ 2, noise(data))
    vis = measurement(data)
    lnorm = VLBILikelihoods.lognorm(CoherencyLikelihood(vis, Σ))
    ℓ = ConditionedLikelihood(_Coherency(Σ, lnorm), vis)
    return ℓ
end

struct _VisAmp{S}
    S::S
end

function (v::_VisAmp)(μ)
    return RiceAmplitudeLikelihood(abs.(baseimage(μ)), v.S)
end

# internal function that creates the likelihood for a set of visibility amplitudes
function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTVisibilityAmplitudeDatum})
    Σ = noise(data) .^ 2
    amp = measurement(data)
    ℓ = ConditionedLikelihood(_VisAmp(Σ), amp)
    return ℓ
end

struct _LCamp{F, S, L}
    f::F
    S::S
    L::L
end

function (c::_LCamp)(μ)
    return AmplitudeLikelihood(c.f(baseimage(μ)), c.S, c.L)
end

# internal function that creates the likelihood for a set of log closure amplitudes
function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTLogClosureAmplitudeDatum})
    Σlca = factornoisecovariance(arrayconfig(data))
    f = Base.Fix2(logclosure_amplitudes, designmat(arrayconfig(data)))
    amp = measurement(data)
    lnorm = VLBILikelihoods.lognorm(AmplitudeLikelihood(amp, Σlca))
    ℓ = ConditionedLikelihood(_LCamp(f, Σlca, lnorm), amp)
    return ℓ
end

struct _CPhase{F, S, L}
    f::F
    S::S
    L::L
end

function (c::_CPhase)(μ)
    return ClosurePhaseLikelihood(c.f(baseimage(μ)), c.S, c.L)
end

# internal function that creates the likelihood for a set of closure phase datum
function makelikelihood(data::Comrade.EHTObservationTable{<:Comrade.EHTClosurePhaseDatum})
    Σcp = factornoisecovariance(arrayconfig(data))
    f = Base.Fix2(closure_phases, designmat(arrayconfig(data)))
    phase = measurement(data)
    lnorm = VLBILikelihoods.lognorm(ClosurePhaseLikelihood(phase, Σcp))
    ℓ = ConditionedLikelihood(_CPhase(f, Σcp, lnorm), phase)
    return ℓ
end

export ImgNormalData

"""
    ImgNormalData(reduction, measurement, noise)

Container for image-domain data assuming a normal likelihood with a `measurement` and `noise`. 
"""
struct ImgNormalData{F, M, N}
    reduction::F
    measurement::M
    noise::N
end


struct _Reduced{F, T}
    f::F
    σ::T
end

struct NormalFast{T, S}
    μ::T
    σ::S
end

function DensityInterface.logdensityof(d::NormalFast, x)
    dev = (x .- d.μ) .* inv.(d.σ)
    exponent = -sum(abs2, dev) / 2
    return exponent
end

function (c::_Reduced)(μ)
    return NormalFast(c.f(μ), c.σ)
end

function makelikelihood(data::ImgNormalData)
    σ = data.noise
    meas = data.measurement
    f = data.reduction
    ℓ = ConditionedLikelihood(_Reduced(f, σ), meas)
    return ℓ
end
