export transform, asflat, ascube

"""
    $(TYPEDEF)
A transformed version of a `VLBIPosterior` object.
This is an internal type that an end user shouldn't have to directly construct.
To construct a transformed posterior see the [`asflat`](@ref asflat), [`ascube`](@ref ascube).
"""
struct TransformedVLBIPosterior{P <: VLBIPosterior, T} <: AbstractVLBIPosterior
    lpost::P
    transform::T
end
(post::TransformedVLBIPosterior)(θ) = logdensityof(post, θ)
admode(post::TransformedVLBIPosterior) = admode(post.lpost)

# Latent-space tags for a `TransportedDistribution`: `TVFlat` (unconstrained ℝⁿ, `asflat`)
# has `stop === nothing`; the unit hypercube (`ascube`) carries a `StdUniform`. These are the
# single place that encodes PT's type-parameter layout — dispatch on them, don't respell them.
const FlatTransport = TransportedDistribution{<:Any, <:Any, Nothing}
const CubeTransport = TransportedDistribution{<:Any, <:Any, <:StdUniform}

# Is the transformed posterior over the unit hypercube ([0,1]^n, i.e. `ascube`)?
_is_cube(t::TransportedDistribution) = t isa CubeTransport
_is_cube(post::TransformedVLBIPosterior) = _is_cube(post.transform)

function prior_sample(rng, tpost::TransformedVLBIPosterior, args...)
    inv = Base.Fix1(latent_pback, tpost)
    return map(inv, prior_sample(rng, tpost.lpost, args...))
end
function prior_sample(rng, tpost::TransformedVLBIPosterior)
    inv = Base.Fix1(latent_pback, tpost)
    return inv(prior_sample(rng, tpost.lpost))
end

dimension(post::TransformedVLBIPosterior) = dimension(post.transform)


_transport_to_post(post::VLBIPosterior, space) =
    TransformedVLBIPosterior(post, transport_to(post.prior, space))

"""
    transport_to(post::VLBIPosterior, space)

Build a [`TransformedVLBIPosterior`](@ref) whose parameters live in the latent `space`
(`TVFlat()` for unconstrained ℝⁿ, `StdUniform()` for the unit hypercube). The
[`asflat`](@ref asflat) / [`ascube`](@ref ascube) helpers call this with the respective
spaces.

Also accepts an already-transformed posterior: it is returned unchanged when it already
lives in `space`, and rebuilt from the underlying `VLBIPosterior` otherwise — an explicit
space request is always honored.

The space argument is typed (rather than left as `::Any`) so these methods stay strictly
more specific than ProbabilityTransports' own `transport_to(dist, ::AbstractStdDist)` /
`transport_to(dist, ::TVFlat)` and don't collide with them by ambiguity.
"""
PT.transport_to(post::VLBIPosterior, space::PT.AbstractStdDist) = _transport_to_post(post, space)
PT.transport_to(post::VLBIPosterior, space::PT.TVFlat) = _transport_to_post(post, space)

# Does the transported distribution already live in `space`? Matched by latent family
# (`StdUniform` vs `StdNormal` vs flat), ignoring eltype/dimension parameters — hence the
# comparison against the space's unparameterized wrapper type.
_same_latent(t::TransportedDistribution, ::PT.TVFlat) = t isa FlatTransport
_same_latent(t::TransportedDistribution, space::PT.AbstractStdDist) =
    getfield(t, :stop) isa Base.typename(typeof(space)).wrapper

_retransport(tpost::TransformedVLBIPosterior, space) =
    _same_latent(tpost.transform, space) ? tpost : transport_to(tpost.lpost, space)
PT.transport_to(tpost::TransformedVLBIPosterior, space::PT.AbstractStdDist) = _retransport(tpost, space)
PT.transport_to(tpost::TransformedVLBIPosterior, space::PT.TVFlat) = _retransport(tpost, space)

"""
    maybe_transport(post, space)

Backend helper for a sampler's `transport_method`-style keyword. `space === nothing` means
"use the backend default": `asflat` for a raw [`VLBIPosterior`](@ref), pass-through for an
already-transformed posterior. An explicit `space` is always honored via
[`transport_to`](@ref transport_to), including re-transporting a `TransformedVLBIPosterior`.
"""
maybe_transport(post::VLBIPosterior, ::Nothing) = asflat(post)
maybe_transport(tpost::TransformedVLBIPosterior, ::Nothing) = tpost
maybe_transport(post::AbstractVLBIPosterior, space) = transport_to(post, space)


"""
    transform(posterior::TransformedVLBIPosterior, x)

Transforms the value `x` from the transformed space (e.g. unit hypercube if using [`ascube`](@ref ascube))
to parameter space which is usually encoded as a `NamedTuple`.

For the inverse transform see [`inverse`](@ref inverse)
"""
PT.latent_pfwd(p::TransformedVLBIPosterior, x) = latent_pfwd(p.transform, x)
transform(p::TransformedVLBIPosterior, x) = latent_pfwd(p, x)


"""
    inverse(posterior::TransformedVLBIPosterior, x)

Transforms the value `x` from parameter space to the transformed space
(e.g. unit hypercube if using [`ascube`](@ref ascube)).

For the forward transform see [`transform`](@ref transform)
"""
PT.latent_pback(p::TransformedVLBIPosterior, x) = latent_pback(p.transform, x)
inverse(p::TransformedVLBIPosterior, x) = latent_pback(p, x)

"""
    asflat(post::VLBIPosterior)

Construct a flattened version of the posterior where the parameters are transformed to live in
(-∞, ∞).

This returns a `TransformedVLBIPosterior` that obeys the `DensityInterface` and can be evaluated
in the usual manner, i.e. `logdensityof`. Note that the transformed posterior automatically
includes the terms log-jacobian terms of the transformation.

# Example
```julia-repl
julia> tpost = asflat(post)
julia> x0 = prior_sample(tpost)
julia> logdensityof(tpost, x0)
```

# Notes
This is the transform that should be used if using typical MCMC methods, i.e. NUTS.
For the transformation to the unit hypercube see [`ascube`](@ref ascube)
"""
asflat(post::VLBIPosterior) = transport_to(post, TVFlat())

"""
    ascube(post::VLBIPosterior)

Construct a flattened version of the posterior where the parameters are transformed to live in
(0, 1), i.e. the unit hypercube.

This returns a `TransformedVLBIPosterior` that obeys the `DensityInterface` and can be evaluated
in the usual manner, i.e. `logdensityof`. Note that the transformed posterior automatically
includes the terms log-jacobian terms of the transformation.

# Example
```julia-repl
julia> tpost = ascube(post)
julia> x0 = prior_sample(tpost)
julia> logdensityof(tpost, x0)
```

# Notes
This is the transform that should be used if using typical NestedSampling methods,
i.e. `ComradeNested`. For the transformation to unconstrained space see [`asflat`](@ref asflat)
"""
ascube(post::VLBIPosterior) = transport_to(post, StdUniform())


function Base.show(io::IO, mime::MIME"text/plain", post::TransformedVLBIPosterior)
    println(io, "TransformedVLBIPosterior(")
    show(io, mime, post.lpost)
    if _is_cube(post)
        println(io, "Transform: Params to [0,1]^$(dimension(post))")
    else
        println(io, "Transform: Params to ℝ^$(dimension(post))")
    end
    return print(io, ")")
end


# A single log-density covers both the flat (`TVFlat`) and cube (`StdUniform`) cases.
# `latent_pfwd_and_logdensity` returns the transported point `p` together with the
# pulled-back *prior* log-density `ℓ`:
#   * flat:  ℓ = logpdf(prior, p) + logjac  ⇒  full posterior with the likelihood;
#   * cube:  ℓ = logpdf(StdUniform, x) = 0 inside [0,1]^n and -Inf outside ⇒ the prior
#            is absorbed by the map, so only the likelihood (bounds-checked) remains.
@inline function DensityInterface.logdensityof(post::TransformedVLBIPosterior, x::AbstractArray)
    p, ℓ = latent_pfwd_and_logdensity(post.transform, x)
    return loglikelihood(post.lpost, p) + ℓ
end

# Cube specialization: check the bounds (`logpdf(StdUniform, x)`, no transport) *before*
# running the transport and forward model, so out-of-cube proposals (e.g. slice-sampler
# expansion steps) short-circuit to -Inf instead of evaluating the full likelihood.
@inline function DensityInterface.logdensityof(
        post::TransformedVLBIPosterior{P, <:CubeTransport}, x::AbstractArray
    ) where {P <: VLBIPosterior}
    ℓ = Dists.logpdf(post.transform, x)
    isfinite(ℓ) || return ℓ
    return loglikelihood(post.lpost, latent_pfwd(post.transform, x)) + ℓ
end
