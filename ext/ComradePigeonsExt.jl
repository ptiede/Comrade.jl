module ComradePigeonsExt

using Comrade

using Pigeons
using AbstractMCMC
using LogDensityProblems
using HypercubeTransform
using TransformVariables
using Random


Pigeons.initialization(tpost::Comrade.TransformedVLBIPosterior, rng::Random.AbstractRNG, ::Int) = prior_sample(rng, tpost)

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

Pigeons.default_explorer(::Comrade.TransformedVLBIPosterior{P,<:HypercubeTransform.AbstractHypercubeTransform}) where {P}  =
    SliceSampler()

Pigeons.default_explorer(::Comrade.TransformedVLBIPosterior{P,<:TransformVariables.AbstractTransform}) where {P} =
    Pigeons.AutoMALA(;default_autodiff_backend = :Zygote)

function Pigeons.default_reference(tpost::Comrade.TransformedVLBIPosterior)
    t = tpost.transform
    p = tpost.lpost.prior
    return PriorRef(p, t)
end

function Pigeons.sample_iid!(target::Comrade.TransformedVLBIPosterior, replica, shared)
    replica.state = Pigeons.initialization(target, replica.rng, replica.replica_index)
end

function Pigeons.sample_iid!(target::PriorRef{P, <:TransformVariables.AbstractTransform}, replica, shared) where {P}
    replica.state .= Comrade.inverse(target.transform, rand(replica.rng, target.model))
end

function Pigeons.sample_iid!(::PriorRef{P,<:HypercubeTransform.AbstractHypercubeTransform}, replica, shared) where {P}
    rand!(replica.rng, replica.state)
end


function Pigeons.sample_array(tpost::Comrade.TransformedVLBIPosterior, pt::Pigeons.PT)
    samples = sample_array(pt)
    tbl = mapreduce(hcat, eachslice(samples, dims=(3,), drop=true)) do arr
        s = map(x->Comrade.transform(tpost, @view(x[begin:end-1])), eachrow(arr))
        return s
    end

    sts = (logdensity= samples[:, end, :] |> vec,)

    return Comrade.PosteriorSamples(tbl, sts; metadata=Dict(:sampler=>:Pigeons, :post=>tpost))
end


LogDensityProblems.dimension(t::PriorRef) = Comrade.dimension(t.transform)
LogDensityProblems.logdensity(t::PriorRef, x) = t(x)
LogDensityProblems.capabilities(::Type{<:PriorRef}) = LogDensityProblems.LogDensityOrder{0}()

end
