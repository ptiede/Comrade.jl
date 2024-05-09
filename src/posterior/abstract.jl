using LinearAlgebra
using VLBILikelihoods

export vlbimodel, logdensityof, dimension, skymodel, instrumentmodel, dataproducts,
       forward_model, prior_sample, simulate_observation,
       VLBIPosterior, logdensityof, loglikelihood

abstract type AbstractVLBIPosterior end
@inline DensityInterface.DensityKind(::AbstractVLBIPosterior) = DensityInterface.IsDensity()
function logprior(d::AbstractVLBIPosterior, θ)
    # @info "HERE"
    logdensityof(d.prior, θ)
end
LogDensityProblems.logdensity(d::AbstractVLBIPosterior, θ) = logdensityof(d, θ)
LogDensityProblems.dimension(d::AbstractVLBIPosterior) = dimension(d)
LogDensityProblems.capabilities(::Type{<:AbstractVLBIPosterior}) = LogDensityProblems.LogDensityOrder{0}()
skymodel(d::AbstractVLBIPosterior) = getfield(d, :skymodel)
instrumentmodel(d::AbstractVLBIPosterior) = getfield(d, :instrumentmodel)
HypercubeTransform.dimension(d::AbstractVLBIPosterior) = length(d.prior)

@noinline logprior_ref(d, x) = logprior(d, x[])

function ChainRulesCore.rrule(::typeof(logprior), d::AbstractVLBIPosterior, x)
    p = logprior(d, x)
    # We need this
    function _logprior_pullback(Δ)
        # @info "HERE"
        xr = Ref(x)
        dxr = Ref(ntzero(x))
        autodiff(Reverse, logprior_ref, Active, Const(d), Duplicated(xr, dxr))
        return NoTangent(), NoTangent(), Comrade.rmap(x->Δ*x, dxr[])
    end
    return p, _logprior_pullback
end

function DensityInterface.logdensityof(post::AbstractVLBIPosterior, x)
    pr = logprior(post, x)
    !isfinite(pr) && return -Inf
    return loglikelihood(post, x) + pr
end

prior_sample(post::AbstractVLBIPosterior, args...) = prior_sample(Random.default_rng(), post, args...)

"""
    prior_sample([rng::AbstractRandom], post::AbstractVLBIPosterior)

Returns a single sample from the prior distribution.
"""
function prior_sample(rng, post::AbstractVLBIPosterior)
    return rand(rng, post.prior)
end


function forward_model(d::AbstractVLBIPosterior, θ)
    vis = idealvisibilities(skymodel(d), θ)
    return apply_instrument(vis, instrumentmodel(d), θ)
end


function loglikelihood(d::AbstractVLBIPosterior, θ)
    vis = forward_model(d, θ)
    # Convert because of conventions
    return logdensityofvis(d.lklhds, vis)
end

"""
    skymodel(post::AbstractVLBIPosterior, θ)

Returns the sky model or image of a posterior using the parameter values`θ`
"""
function skymodel(post::AbstractVLBIPosterior, θ)
    return skymodel(post.skymodel, θ.sky)
end

"""
    instrumentmodel(post::AbstractVLBIPosterior, θ; decompose=true)

Returns the instrument model for the posterior with parameters `θ`.
If `decompose` is true, then each Jones matrix of the instrument model is returned individually.
Note that for Stokes I imaging the complex gains are returned.
"""
function instrumentmodel(post::AbstractVLBIPosterior, θ; decompose=true)
    jonesmatrices(instrumentmodel(post), θ; decompose)
end

"""
    dataproducts(d::AbstractRadioLikelihood)

Returns the data products you are fitting as a tuple. The order of the tuple corresponds
to the order of the `dataproducts` argument in [`AbstractRadioLikelihood`](@ref).
"""
function dataproducts(d::AbstractVLBIPosterior)
    return getfield(d, :data)
end



function logdensityofvis(lklhds, vis::AbstractArray)
    fl = Base.Fix2(logdensityof, vis)
    ls = map(fl, lklhds)
    return sum(ls)
end


include("likelihood.jl")
include("vlbiposterior.jl")
include("transformed.jl")
