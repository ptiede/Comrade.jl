using LinearAlgebra
using VLBILikelihoods: VLBILikelihoods, AmplitudeLikelihood, ClosurePhaseLikelihood,
    CoherencyLikelihood, ComplexVisLikelihood, RiceAmplitudeLikelihood
# One `loglikelihood` with Distributions and StatsBase, so loading them with Comrade does not clash.
import StatsBase: loglikelihood

export logdensityof, dimension, skymodel, instrumentmodel, dataproducts,
    forward_model, prior_sample, simulate_observation,
    VLBIPosterior, loglikelihood, chi2, residuals, ndata


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
logprior(d::AbstractVLBIPosterior, ::NamedTuple{()}) = zero(real(eltype(baseimage(weight(_commonvis(d.data))))))

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
dimension(d::AbstractVLBIPosterior) = length(d.prior)
EnzymeRules.inactive(::typeof(instrumentmodel), args...) = nothing


function DensityInterface.logdensityof(post::AbstractVLBIPosterior, x)
    pr = logprior(post, x)
    l = loglikelihood(post, x)
    # isnan(l) && throw(ArgumentError("NaN in loglikelihood at $x"))
    return l + pr
end

prior_sample(post::AbstractVLBIPosterior, args...) = prior_sample(Random.default_rng(), post, args...)

"""
    prior_sample([rng::AbstractRNG], post::AbstractVLBIPosterior, dims...)

Returns a sample from the prior distribution of the posterior. If `dims` are given, returns an
array of size `dims` of independent draws.
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
    img, vis = idealmaps(D, skymodel(d), θ)
    return img, apply_instrument(vis, instrumentmodel(d), _instrumentvalues(θ))
end

# The instrument parameters of `θ`; an ideal instrument has none.
_instrumentvalues(θ::NamedTuple) = haskey(θ, :instrument) ? θ.instrument : (;)

@inline function forward_model(d::AbstractVLBIPosterior, θ)
    return forward_model_map(datatype(typeof(d)), d, θ)
end

"""
    loglikelihood(d::AbstractVLBIPosterior, θ)

Computes the log-likelihood of the posterior `d` with parameters `θ`.
"""
@inline function loglikelihood(d::AbstractVLBIPosterior, θ)
    img, vis = forward_model(d, θ)
    # Convert because of conventions
    lis = d.lklhdsimg
    return logdensityofvis(d.lklhds, vis) + logdensityofimg(lis, img, θ)
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

The Jones matrices `J_a U_a` of a posterior's instrument at the parameter values `θ`, for every
station, time stamp and channel of its data, as an array over `(Ant, Ti, Fr)` (see
[`forward_jones`](@ref)). A multi-field posterior gives one array per field, keyed by field name.
"""
function instrumentmodel(post::AbstractVLBIPosterior, θ)
    return _forwardjones(instrumentmodel(post), θ, post.data)
end

_forwardjones(J, θ, data::Tuple) = forward_jones(_observedmodel(J), _instrumentvalues(θ), _commonvis(data))
_forwardjones(J::FieldInstruments, θ, data::NamedTuple) = map((j, d) -> _forwardjones(j, θ, d), J.fields, data)

_observedmodel(J::ObservedInstrumentModel) = J
_observedmodel(J::ObservedSelfCalInstrumentModel) = J.model


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

@inline _imglogdensity(l, img, θ) = logdensityof(l, img)

@inline function logdensityofimg(lklhds::Tuple, img, θ)
    ls = map(l -> _imglogdensity(l, img, θ), lklhds)
    return sum(ls)
end

@inline logdensityofimg(lklhds::Tuple{}, img, θ) = logdensityofimg(lklhds, img)

## There is no image data so just return 0
@inline function logdensityofimg(lklhds::Tuple{}, img::Number)
    return zero(img)
end

@inline function logdensityofimg(lklhds::Tuple{}, img::AbstractArray{T}) where {T}
    return zero(T)
end

# The components of a `MultiSkyModel`.
@inline logdensityofimg(lklhds::Tuple{}, img::Tuple) = sum(map(i -> logdensityofimg(lklhds, i), img))


include("likelihood.jl")
include("vlbiposterior.jl")
include("fields.jl")
include("transformed.jl")
