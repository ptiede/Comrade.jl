export transform, asflat, ascube

"""
    $(TYPEDEF)
A transformed version of a `VLBIPosterior` object.
This is an internal type that an end user shouldn't have to directly construct.
To construct a transformed posterior see the [`asflat`](@ref asflat), [`ascube`](@ref ascube).
"""
struct TransformedVLBIPosterior{P<:VLBIPosterior,T} <: AbstractVLBIPosterior
    lpost::P
    transform::T
end
(post::TransformedVLBIPosterior)(θ) = logdensityof(post, θ)


function prior_sample(rng, tpost::TransformedVLBIPosterior, args...)
    inv = Base.Fix1(HypercubeTransform.inverse, tpost)
    map(inv, prior_sample(rng, tpost.lpost, args...))
end
function prior_sample(rng, tpost::TransformedVLBIPosterior)
    inv = Base.Fix1(HypercubeTransform.inverse, tpost)
    inv(prior_sample(rng, tpost.lpost))
end

HypercubeTransform.dimension(post::TransformedVLBIPosterior) = dimension(post.transform)


"""
    transform(posterior::TransformedVLBIPosterior, x)

Transforms the value `x` from the transformed space (e.g. unit hypercube if using [`ascube`](@ref ascube))
to parameter space which is usually encoded as a `NamedTuple`.

For the inverse transform see [`inverse`](@ref HypercubeTransform.inverse)
"""
HypercubeTransform.transform(p::TransformedVLBIPosterior, x) = transform(p.transform, x)


"""
    inverse(posterior::TransformedVLBIPosterior, x)

Transforms the value `y` from parameter space to the transformed space
(e.g. unit hypercube if using [`ascube`](@ref ascube)).

For the inverse transform see [`transform`](@ref HypercubeTransform.transform)
"""
HypercubeTransform.inverse(p::TransformedVLBIPosterior, y) = HypercubeTransform.inverse(p.transform, y)

"""
    asflat(post::VLBIPosterior)

Construct a flattened version of the posterior where the parameters are transformed to live in
(-∞, ∞).

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
This is the transform that should be used if using typical MCMC methods, i.e. `ComradeAHMC`.
For the transformation to the unit hypercube see [`ascube`](@ref ascube)

"""
function HypercubeTransform.asflat(post::VLBIPosterior)
    pr = post.prior
    tr = asflat(pr)
    return TransformedVLBIPosterior(post, tr)
end

function Base.show(io::IO, mime::MIME"text/plain", post::TransformedVLBIPosterior{P, T}) where {P, T<:TV.AbstractTransform}
    println(io, "TransformedVLBIPosterior(")
    show(io, mime, post.lpost)
    println(io, "Transform: Params to ℝ^$(dimension(post))")
    print(io, ")")
end

# @noinline function _refed_transform_and_logprior(t::TV.AbstractTransform, prior, x)
#     # @info "HERE"
#     p, logjac = transform_and_logjac(t, x)
#     logpr = logdensityof(prior, p)
#     return p, logpr + logjac
# end

# @noinline function _refered_transform_and_logprior!(y, t, prior, x)
#     out = _refed_transform_and_logprior(t, prior, x)
#     y[] = out
#     nothing
# end

# function ChainRulesCore.rrule(::typeof(_refed_transform_and_logprior), t::TV.AbstractTransform, prior, x)
#     out = _refed_transform_and_logprior(t, prior, x)
#     px = ProjectTo(x)
#     function _refed_transform_and_logprior_pullback(Δ)
#         # @info typeof(Δ)
#         Δout = make_sensitivity(out, Δ)
#         @info Δout
#         # @info typeof(out)
#         Δoutr = Ref(Δout)
#         outr = Ref(ntzero(out))
#         dx = zero(x)
#         autodiff(Reverse, _refered_transform_and_logprior!, Duplicated(outr, Δoutr), Const(t), Const(prior), Duplicated(x, dx))
#         return NoTangent(), NoTangent(), NoTangent(), px(dx)
#     end
#     return out, _refed_transform_and_logprior_pullback
# end

# function make_sensitivity(out::NTuple{N, Any}, Δ::Tangent) where {N}
#     # @info typeof(Δ)
#     return ntuple(i->make_sensitivity(out[i], Δ[i]), Val(N))
# end

# function make_sensitivity(out::NamedTuple{N}, Δ::Tangent) where {N}
#     # @info typeof(Δ)
#     return NamedTuple{N}(map(N) do N
#         return make_sensitivity(out[N], Δ[N])
#     end)
# end

# make_sensitivity(out, Δ::Tangent) = throw(ArgumentError("wtf $(typeof(Δ))"))
# make_sensitivity(out, Δ) = (Δ)
# make_sensitivity(out, ::AbstractZero) = zero(out)


@inline function DensityInterface.logdensityof(post::TransformedVLBIPosterior{P, T}, x::AbstractArray) where {P, T<:TV.AbstractTransform}
    p, logjac = transform_and_logjac(post.transform, x)
    pr = logprior(post.lpost, p)
    lp = post.lpost
    return loglikelihood(lp, p) + pr + logjac
end


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
function HypercubeTransform.ascube(post::VLBIPosterior)
    pr = post.prior
    tr = ascube(pr)
    return TransformedVLBIPosterior(post, tr)
end

function Base.show(io::IO, post::TransformedVLBIPosterior{P, T}) where {P, T<:HypercubeTransform.AbstractHypercubeTransform}
    println(io, "TransformedVLBIPosterior(")
    println(io, post.lpost)
    println(io, "Transform: Params to [0,1]^$(dimension(post))")
    print(io, ")")
end

function DensityInterface.logdensityof(tpost::TransformedVLBIPosterior{P, T}, x::AbstractArray) where {P, T<:HypercubeTransform.AbstractHypercubeTransform}
    # Check that x really is in the unit hypercube. If not return -Inf
    for xx in x
        (xx > 1 || xx < 0) && return -Inf
    end
    p = transform(tpost.transform, x)
    post = tpost.lpost
    return loglikelihood(post, p)
end
