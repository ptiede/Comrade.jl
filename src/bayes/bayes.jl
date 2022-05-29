export Posterior, prior_sample, asflat, ascube, flatten, logdensityof, transform, inverse, dimension

import DensityInterface
import ParameterHandling
import ParameterHandling: flatten
using HypercubeTransform
using TransformVariables
using ValueShapes: NamedTupleDist

"""
    Posterior(lklhd, prior, model)
Posterior density that follows obeys the [DensityInferface](https://github.com/JuliaMath/DensityInferface.jl)

The expected arguments are:

- lklhd: Which should be an intance of RadioLikelihood with whatever data products you want to fit
- prior: This should be a `NamedTuple` with the priors for each model ParameterHandling
- model: Function that takes a `NamedTuple` of parameters and constructs the `Comrade` model.

To evaluate the logdensity of the posterior use can either use `DensityInterface.logdensityof`
or `logdensity`.
"""
struct Posterior{L,P,F}
    lklhd::L
    prior::P
    model::F
end

function Posterior(lklhd, prior::NamedTuple, model)
    return Posterior(lklhd, NamedTupleDist(prior), model)
end

@inline DensityInterface.DensityKind(::Posterior) = IsDensity()

function DensityInterface.logdensityof(post::Posterior, x)
    pr = logdensityof(post.prior, x)
    !isfinite(pr) && return -Inf
    return logdensityof(post.lklhd, post.model(x)) + pr
end

"""
    $(SIGNATURES)
Samples the prior distribution from the posterior `nsamples` times.

Returns a Vector{NamedTuple} that can be used to initialize optimization and sample algorithms
"""
function prior_sample(post::Posterior, args...)
    return rand(post.prior, args...)
end

function prior_sample(post::Posterior)
    return rand(post.prior)
end

"""
    `$(TYPEDEF)`
A transformed version of the posterior `lpost`. This is an internal type
that an end user shouldn't have to directly construct.

To construct a transformed posterior see the `asflat`, `ascube`, and `flatten`
functions.

When evaluating the `logdensity` of the `TransformedPosterior`, `logdenisty`
will include any jacobian terms automatically, and expects the argument to be the
vector of transformed quantities. This is usually easier to incorporate with
samplers, and optimizers, which expect homogenous vectors as arguments.
"""
struct TransformedPosterior{P<:Posterior,T}
    lpost::P
    transform::T
end

function prior_sample(tpost::TransformedPosterior, args...)
    inv = Base.Fix1(HypercubeTransform.inverse, post)
    map(inv, prior_sample(tpost.lpost, args...))
end

function prior_sample(tpost::TransformedPosterior)
    inv = Base.Fix1(HypercubeTransform.inverse, post)
    inv(prior_sample(tpost.lpost))
end



"""
    transform(posterior::TransformedPosterior, x)
Transforms the value `x` into parameter space, i.e. usually a NamedTuple.
"""
HypercubeTransform.transform(p::TransformedPosterior, x) = transform(p.transform, x)


"""
    inverse(posterior::TransformedPosterior, x)
Transforms the value model parameters `x` into the flattened transformed space.
"""
HypercubeTransform.inverse(p::TransformedPosterior, y) = HypercubeTransform.inverse(p.transform, y)


@inline DensityInterface.DensityKind(::TransformedPosterior) = DensityInterface.IsDensity()

# MeasureBase.logdensityof(tpost::Union{Posterior,TransformedPosterior}, x) = DensityInterface.logdensityof(tpost, x)


"""
    asflat(post::Posterior)
Construct a flattened version of the posterior, where the parameters are transformed so that
their support is from (-∞, ∞). This uses [TransformVariables](https://github.com/tpapp/TransformVariables.jl)

The transformed posterior can then be evaluated by the `logdensityof(transformed_posterior,x)`
method following the `DensityInterface`, where `x` is a flattened vector of the infinite support
variables. **Note** this already includes the jacobian of the transformation so this does
not need to be added.

This is useful for optimization and sampling algorithms such as HMC that will use gradients
to explore the posterior surface.
"""
function HypercubeTransform.asflat(post::Posterior)
    pr = getfield(post.prior, :_internal_distributions)
    tr = asflat(pr)
    return TransformedPosterior(post, tr)
end

@inline function DensityInterface.logdensityof(post::TransformedPosterior{P, T}, x::AbstractArray) where {P, T<:HypercubeTransform.NamedFlatTransform}
    p, logjac = transform_and_logjac(post.transform, x)
    return logdensityof(post.lpost, p) + logjac
end

HypercubeTransform.dimension(post::TransformedPosterior) = dimension(post.transform)
HypercubeTransform.dimension(post::Posterior) = dimension(asflat(post))

"""
    ascube(post::Posterior)
Construct a flattened version of the posterior, where the parameters are transformed to live
in the unit hypercube. In astronomy parlance, we are transforming the variables to the unit
hypercube. This is done using the [HypercubeTransform](https://github.com/ptiede/HypercubeTransform.jl)
package.

The transformed posterior can then be evaluated by the `logdensityof(transformed_posterior,x)`
method following the `DensityInterface`, where `x` vector that lives in the unit hypercube.
**Note** this already includes the jacobian of the transformation so this does
not need to be added.

This transform is useful for NestedSampling methods that often assume that the model is written
to live in the unit hypercube.
"""
function HypercubeTransform.ascube(post::Posterior)
    pr = getfield(post.prior, :_internal_distributions)
    tr = ascube(pr)
    return TransformedPosterior(post, tr)
end

function DensityInterface.logdensityof(tpost::TransformedPosterior{P, T}, x::AbstractArray) where {P, T<:HypercubeTransform.AbstractHypercubeTransform}
    # Check that x really is in the unit hypercube. If not return -Inf
    for xx in x
        (xx > 1 || xx < 0) && return -Inf
    end
    p = transform(tpost.transform, x)
    post = tpost.lpost
    return logdensityof(post.lklhd, post.model(p))
end

struct FlatTransform{T}
    transform::T
end

HypercubeTransform.transform(t::FlatTransform, x) = t.transform(x)
HypercubeTransform.inverse(::FlatTransform, x) = first(ParameterHandling.flatten(x))
HypercubeTransform.dimension(t::FlatTransform) = t.transform.unflatten.sz[end]
"""
    flatten(post::Posterior)
Flatten the representation of the posterior. Internally this uses ParameterHandling to
construct a flattened version of the posterior.

Note this is distinct from `asflat` that transforms the variables to live in (-∞,∞).
Instead this method just flattens the repsentation of the model from a NamedTuple to a vector.
This allows the easier integration to optimization and sampling algorithms.
"""
function ParameterHandling.flatten(post::Posterior)
    x0 = rand(post.prior)
    _, unflatten = ParameterHandling.flatten(x0)
    return TransformedPosterior(post, FlatTransform(unflatten))
end

# function DensityInterface.logdensityof(post::TransformedPosterior{P,T}, x) where {P, T<: HypercubeTransform.NamedFlatTransform}
#     return logdensity(post.lpost, transform(post.transform, x))
# end
