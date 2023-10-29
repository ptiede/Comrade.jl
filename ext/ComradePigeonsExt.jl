module ComradePigeonsExt

using Comrade

if isdefined(Base, :get_extension)
    using Pigeons
    using AbstractMCMC
    using LogDensityProblems
    using HypercubeTransform
    using TypedTables
    using TransformVariables
    using Random

else
    using ..Pigeons
    using ..AbstractMCMC
    using ..LogDensityProblems
    using ..HypercubeTransform
    using ..TypedTables
    using ..TransformVariables
    using ..Random
end

Pigeons.initialization(tpost::Comrade.TransformedPosterior, rng::Random.AbstractRNG, ::Int) = prior_sample(rng, tpost)

struct PriorRef{P,T}
    model::P
    transform::T
end

function (p::PriorRef{P,<:TransformVariables.AbstractTransform})(x) where {P}
    y, lj = TransformVariables.transform_and_logjac(p.transform, x)
    logdensityof(p.model, y) + lj
end

function (p::PriorRef{P,<:HypercubeTransform.AbstractHypercubeTransform})(x) where {P}
    for xx in x
        (xx > 1 || xx < 0) && return convert(eltype(x), -Inf)
    end
    return zero(eltype(x))
end

Pigeons.default_explorer(::Comrade.TransformedPosterior{P,<:HypercubeTransform.AbstractHypercubeTransform}) where {P}  =
    SliceSampler()

Pigeons.default_explorer(::Comrade.TransformedPosterior{P,<:TransformVariables.AbstractTransform}) where {P} =
    Pigeons.AutoMALA(;default_autodiff_backend = :Zygote)

function Pigeons.default_reference(tpost::Comrade.TransformedPosterior)
    t = tpost.transform
    p = tpost.lpost.prior
    return PriorRef(p, t)
end

function Pigeons.sample_iid!(target::Comrade.TransformedPosterior, replica, shared)
    replica.state = Pigeons.initialization(target, replica.rng, replica.replica_index)
end

function Pigeons.sample_iid!(target::PriorRef{P, <:TransformVariables.AbstractTransform}, replica, shared) where {P}
    replica.state .= Comrade.inverse(target.transform, rand(replica.rng, target.model))
end

function Pigeons.sample_iid!(::PriorRef{P,<:HypercubeTransform.AbstractHypercubeTransform}, replica, shared) where {P}
    rand!(replica.rng, replica.state)
end


function Pigeons.sample_array(tpost::Comrade.TransformedPosterior, pt::Pigeons.PT)
    samples = sample_array(pt)
    arr =  reshape(samples, size(samples, 1), size(samples, 2))
    return Table(map(x->Comrade.transform(tpost, x), eachrow(arr)))
end


LogDensityProblems.dimension(t::PriorRef) = Comrade.dimension(t.transform)
LogDensityProblems.logdensity(t::PriorRef, x) = t(x)
LogDensityProblems.capabilities(::Type{<:PriorRef}) = LogDensityProblems.LogDensityOrder{0}()

end
