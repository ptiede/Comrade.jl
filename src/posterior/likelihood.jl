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


struct _Visibility{W, L}
    w::W
    L::L
end

(v::_Visibility)(μ) = ComplexVisLikelihood(baseimage(μ), v.w, v.L)

function _vislikelihood(vis, w)
    lnorm = VLBILikelihoods.lognorm(ComplexVisLikelihood(vis, w))
    return ConditionedLikelihood(_Visibility(w, lnorm), vis)
end

struct _Coherency{W, L}
    w::W
    L::L
end

(c::_Coherency)(μ) = CoherencyLikelihood(_storage(baseimage(μ)), c.w, c.L)


function _coherencylikelihood(vis, w)
    lnorm = VLBILikelihoods.lognorm(CoherencyLikelihood(vis, w))
    return ConditionedLikelihood(_Coherency(w, lnorm), vis)
end

struct _VisAmp{W, L}
    w::W
    L::L
end

(v::_VisAmp)(μ) = RiceAmplitudeLikelihood(abs.(baseimage(μ)), v.w, v.L)

function _amplikelihood(amp, w)
    lnorm = VLBILikelihoods.lognorm(RiceAmplitudeLikelihood(amp, w))
    return ConditionedLikelihood(_VisAmp(w, lnorm), amp)
end

struct _LCamp{F, S, L}
    f::F
    S::S
    L::L
end

function (c::_LCamp)(μ)
    return AmplitudeLikelihood(c.f(baseimage(μ)), c.S, c.L)
end

struct _CPhase{F, S, L}
    f::F
    S::S
    L::L
end

function (c::_CPhase)(μ)
    return ClosurePhaseLikelihood(c.f(baseimage(μ)), c.S, c.L)
end

# Coherencies go to the likelihood as dense `(…, 2, 2)` arrays.
function makelikelihood(t::VisTable)
    return _tablelikelihood(eltype(measurement(t)), _dense(measurement(t)), _dense(weight(t)))
end

_tablelikelihood(::Type{<:Complex}, vis, w) = _vislikelihood(vis, w)
_tablelikelihood(::Type{<:Real}, amp, w) = _amplikelihood(amp, w)
_tablelikelihood(::Type{<:StaticArraysCore.StaticMatrix}, vis, w) = _coherencylikelihood(vis, w)

# The closures of model visibilities `vis` over a table's records: their measured values, by
# record then channel, through the design matrix.
struct _ClosureOf{K, D, P}
    keep::K
    design::D
    product::P
end
(c::_ClosureOf)(vis) = c.design * _closurevalue.(Ref(c.product), vis[c.keep])
EnzymeRules.inactive_type(::Type{<:_ClosureOf}) = true

function makelikelihood(ct::ClosureTable)
    Σ = VLBILikelihoods.CholeskyFactor(getfield(ct, :covariance))
    f = _ClosureOf(_measuredvalues(getfield(ct, :vis)), getfield(ct, :design), ct.product)
    x = Vector(baseimage(measurement(ct)))
    if ct.product isa ClosurePhases
        lnorm = VLBILikelihoods.lognorm(ClosurePhaseLikelihood(x, Σ))
        return ConditionedLikelihood(_CPhase(f, Σ, lnorm), x)
    else
        lnorm = VLBILikelihoods.lognorm(AmplitudeLikelihood(x, Σ))
        return ConditionedLikelihood(_LCamp(f, Σ, lnorm), x)
    end
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

export ParamLogDensity

"""
    ParamLogDensity(f)

A log-density term `f(θ)` of the posterior parameters `θ` (the named tuple with `sky` and
`instrument` entries), added to the log-likelihood. Pass it with any image-domain data through
the `imgdata` keyword of [`VLBIPosterior`](@ref).

The term is part of the log-likelihood, not the prior: a transport of the prior to a standard
space does not see it.
"""
struct ParamLogDensity{F}
    f::F
end

makelikelihood(data::ParamLogDensity) = data
@inline _imglogdensity(l::ParamLogDensity, img, θ) = l.f(θ)

function makelikelihood(data::ImgNormalData)
    σ = data.noise
    meas = data.measurement
    f = data.reduction
    ℓ = ConditionedLikelihood(_Reduced(f, σ), meas)
    return ℓ
end
