using LinearAlgebra
using VLBILikelihoods

export vlbimodel, logdensityof, dimension, skymodel, instrumentmodel, dataproducts,
       forward_model, prior_sample, simulate_observation,
       VLBIPosterior, logdensityof, loglikelihood

abstract type AbstractVLBIPosterior end
@inline DensityInterface.DensityKind(::AbstractVLBIPosterior) = DensityInterface.IsDensity()
logprior(d::AbstractVLBIPosterior, θ) = logdensityof(d.prior, θ)
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
    px = ProjectTo(x)
    function _logprior_pullback(Δ)
        # @info "HERE"
        xr = Ref(x)
        dxr = Ref(ntzero(x))
        autodiff(Reverse, logprior_ref, Active, Const(d), Duplicated(xr, dxr))
        return NoTangent(), NoTangent(), (_perturb(Δ, dxr[]))
    end
    return p, _logprior_pullback
end

function _perturb(Δ, x::Union{NamedTuple, Tuple})
    return map(x->_perturb(Δ, x), x)
end

function _perturb(Δ, x)
    return Δ*x
end

function _perturb(Δ, x::AbstractArray)
    x .= Δ*x
    return x
end



function DensityInterface.logdensityof(post::AbstractVLBIPosterior, x)
    pr = logprior(post, x)
    !isfinite(pr) && return -Inf
    return loglikelihood(post, x) + pr
end

prior_sample(post::AbstractVLBIPosterior, args...) = prior_sample(Random.default_rng(), post, args...)

"""
    prior_sample([rng::AbstractRandom], post::AbstractVLBIPosterior, [dims=1])

Returns sample from the prior distribution of the posterior. If dims is specified then
multiple independent draws are returned with shape dims.
"""
function prior_sample(rng, post::AbstractVLBIPosterior)
    return rand(rng, post.prior)
end


function prior_sample(rng, post::AbstractVLBIPosterior, dims)
    map(CartesianIndices(dims)) do _
        return prior_sample(rng, post)
    end
end

"""
    forward_model(d::AbstractVLBIPosterior, θ)

Computes the forward model visibilities of the posterior `d` with parameters `θ`.
Note these are the complex visiblities or the coherency matrices, not the actual
data products observed.
"""
function forward_model(d::AbstractVLBIPosterior, θ)
    vis = idealvisibilities(skymodel(d), θ)
    return apply_instrument(vis, instrumentmodel(d), θ)
end

"""
    loglikelihood(d::AbstractVLBIPosterior, θ)

Computes the log-likelihood of the posterior `d` with parameters `θ`.
"""
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
