using LinearAlgebra
using VLBILikelihoods

export vlbimodel, logdensityof, dimension, skymodel, instrumentmodel, dataproducts,
    forward_model, prior_sample, simulate_observation,
    VLBIPosterior, logdensityof, loglikelihood, chi2, residuals


"""
    $(TYPEDEF)

An abstract VLBI posterior. See [`VLBIPosterior`](@ref) for a concrete implementation.
This implements the  `DensityInterface` and `LogDensityProblem` interfaces.

Default methods include:

 - `logprior(d::AbstractVLBIPosterior, θ)`: Computes the log-prior of the posterior.
 - `loglikelihood(d::AbstractVLBIPosterior, θ)`: Computes the log-likelihood of the posterior.
 - `dimension(d::AbstractVLBIPosterior)`: Returns the dimension of the posterior.
 - `skymodel(d::AbstractVLBIPosterior, θ)`: Returns the sky model of the posterior.
 - `prior_sample(rng::AbstractRandom, d::AbstractVLBIPosterior, dims...)`: Samples from the prior of the posterior.
 - `forward_model(d::AbstractVLBIPosterior, θ)`: Computes the forward model visibilities of the posterior.
 - `dataproducts(d::AbstractVLBIPosterior)`: Returns the data products you are fitting as a tuple.

"""
abstract type AbstractVLBIPosterior end
@inline DensityInterface.DensityKind(::AbstractVLBIPosterior) = DensityInterface.IsDensity()


"""
    logprior(d::AbstractVLBIPosterior, θ)

Computes the log-prior of the posterior `d` with parameters `θ`.
"""
logprior(d::AbstractVLBIPosterior, θ) = logdensityof(d.prior, θ)
logprior(d::AbstractVLBIPosterior, ::NamedTuple{()}) = zero(eltype(d.data[1].config.datatable.U))

LogDensityProblems.logdensity(d::AbstractVLBIPosterior, θ) = logdensityof(d, θ)
LogDensityProblems.dimension(d::AbstractVLBIPosterior) = dimension(d)
LogDensityProblems.capabilities(::Type{<:AbstractVLBIPosterior}) = LogDensityProblems.LogDensityOrder{0}()

"""
    skymodel(d::AbstractVLBIPosterior)

Returns the sky model of the posterior `d`.
"""
skymodel(d::AbstractVLBIPosterior) = getfield(d, :skymodel)

"""
    instrumentmodel(d::AbstractVLBIPosterior)

Returns the instrument model of the posterior `d`.
"""
instrumentmodel(d::AbstractVLBIPosterior) = getfield(d, :instrumentmodel)
HypercubeTransform.dimension(d::AbstractVLBIPosterior) = length(d.prior)
EnzymeRules.inactive(::typeof(instrumentmodel), args...) = nothing


function DensityInterface.logdensityof(post::AbstractVLBIPosterior, x)
    pr = logprior(post, x)
    !isfinite(pr) && return -Inf
    l = loglikelihood(post, x)
    # isnan(l) && throw(ArgumentError("NaN in loglikelihood at $x"))
    return l + pr
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


function prior_sample(rng, post::AbstractVLBIPosterior, dims...)
    return map(CartesianIndices(dims)) do _
        return prior_sample(rng, post)
    end
end

"""
    forward_model(d::AbstractVLBIPosterior, θ)

Computes the forward model visibilities of the posterior `d` with parameters `θ`.
Note these are the complex visiblities or the coherency matrices, not the actual
data products observed.
"""
@inline function forward_model_map(D, d::AbstractVLBIPosterior, θ)
    vis, img = idealmaps(D, skymodel(d), θ)
    return apply_instrument(vis, instrumentmodel(d), θ), img
end

@inline function forward_model(d::AbstractVLBIPosterior, θ)
    forward_model_map(datatype(typeof(d)), d, θ)
end

"""
    loglikelihood(d::AbstractVLBIPosterior, θ)

Computes the log-likelihood of the posterior `d` with parameters `θ`.
"""
@inline function loglikelihood(d::AbstractVLBIPosterior, θ)
    vis, img = forward_model(d, θ)
    # Convert because of conventions
    lis = d.lklhdsimg
    return logdensityofvis(d.lklhds, vis) + logdensityofimg(lis, img)
end

"""
    skymodel(post::AbstractVLBIPosterior, θ)

Returns the sky model or image of a posterior using the parameter values`θ`
"""
function skymodel(post::AbstractVLBIPosterior, θ)
    return skymodel(post.skymodel, θ.sky)
end

"""
    instrumentmodel(post::AbstractVLBIPosterior, θ)

Returns the instrument model of a posterior using the parameter values `θ`.
The output will be a `SiteArray` of the Jones matrices for each site, time, and frequency.
"""
function instrumentmodel(post::AbstractVLBIPosterior, θ)
    J = instrumentmodel(post)
    return forward_jones(J.instrument, θ.instrument)
end


"""
    dataproducts(d::AbstractVLBIPosterior)

Returns the data products you are fitting as a tuple. The order of the tuple corresponds
to the order of the `dataproducts` argument in [`VLBIPosterior`](@ref).
"""
function dataproducts(d::AbstractVLBIPosterior)
    return getfield(d, :data)
end


@inline function logdensityofvis(lklhds, vis::AbstractArray)
    fl = Base.Fix2(logdensityof, vis)
    ls = map(fl, lklhds)
    return sum(ls)
end

@inline function logdensityofimg(lklhds, img::IntensityMap)
    fl = Base.Fix2(logdensityof, img)
    ls = map(fl, lklhds)
    return sum(ls)
end

## There is no image data so just return 0
@inline function logdensityofimg(lklhds, img::Number)
    return zero(img)
end




include("likelihood.jl")
include("vlbiposterior.jl")
include("transformed.jl")
