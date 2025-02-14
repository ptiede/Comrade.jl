module StokedPigeonsExt

using Stoked

using Pigeons
using AbstractMCMC
using LogDensityProblems
using HypercubeTransform
using TransformVariables
using Random


Pigeons.initialization(tpost::Stoked.TransformedVLBIPosterior, rng::Random.AbstractRNG, ::Int) = prior_sample(rng, tpost)

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

Pigeons.default_explorer(::Stoked.TransformedVLBIPosterior{P,<:HypercubeTransform.AbstractHypercubeTransform}) where {P}  =
    SliceSampler()

Pigeons.default_explorer(::Stoked.TransformedVLBIPosterior{P,<:TransformVariables.AbstractTransform}) where {P} =
    Pigeons.AutoMALA(;default_autodiff_backend = :Enzyme)

function Pigeons.default_reference(tpost::Stoked.TransformedVLBIPosterior)
    t = tpost.transform
    p = tpost.lpost.prior
    return PriorRef(p, t)
end

function Pigeons.sample_iid!(target::Stoked.TransformedVLBIPosterior, replica, shared)
    replica.state = Pigeons.initialization(target, replica.rng, replica.replica_index)
end

function Pigeons.sample_iid!(target::PriorRef{P, <:TransformVariables.AbstractTransform}, replica, shared) where {P}
    replica.state .= Stoked.inverse(target.transform, rand(replica.rng, target.model))
end

function Pigeons.sample_iid!(::PriorRef{P,<:HypercubeTransform.AbstractHypercubeTransform}, replica, shared) where {P}
    rand!(replica.rng, replica.state)
end


function Pigeons.sample_array(tpost::Stoked.TransformedVLBIPosterior, pt::Pigeons.PT)
    samples = sample_array(pt)
    tbl = mapreduce(hcat, eachslice(samples, dims=(3,), drop=true)) do arr
        s = map(x->Stoked.transform(tpost, @view(x[begin:end-1])), eachrow(arr))
        return s
    end

    sts = (logdensity= samples[:, end, :] |> vec,)

    return Stoked.PosteriorSamples(tbl, sts; metadata=Dict(:sampler=>:Pigeons, :post=>tpost))
end


LogDensityProblems.dimension(t::PriorRef) = Stoked.dimension(t.transform)
LogDensityProblems.logdensity(t::PriorRef, x) = t(x)
LogDensityProblems.capabilities(::Type{<:PriorRef}) = LogDensityProblems.LogDensityOrder{0}()

end
