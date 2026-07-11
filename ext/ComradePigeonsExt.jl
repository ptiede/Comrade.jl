module ComradePigeonsExt

using Comrade

using Pigeons
using AbstractMCMC
using EnzymeCore
using LogDensityProblems
using Random

const CubeTransport = Comrade.CubeTransport


Pigeons.initialization(tpost::Comrade.TransformedVLBIPosterior, rng::Random.AbstractRNG, ::Int) = prior_sample(rng, tpost)

struct PriorRef{P, T}
    model::P
    transform::T
end

# The Pigeons reference is the prior: `logpdf(::TransportedDistribution, x)` is the
# pulled-back prior log-density in the flat space and the (transport-free) `StdUniform`
# reference density in the cube space.
(p::PriorRef)(x) = Comrade.Dists.logpdf(p.transform, x)

# The cube ([0,1]^n) reference is bounded, so use a slice sampler; every other latent space
# (flat ℝⁿ, StdNormal, …) is unbounded and differentiable, so use gradient-based AutoMALA.
Pigeons.default_explorer(::Comrade.TransformedVLBIPosterior{P, <:CubeTransport}) where {P} =
    SliceSampler()

Pigeons.default_explorer(::Comrade.TransformedVLBIPosterior) =
    Pigeons.AutoMALA(; default_autodiff_backend = :Enzyme)

function Pigeons.default_reference(tpost::Comrade.TransformedVLBIPosterior)
    t = tpost.transform
    p = tpost.lpost.prior
    return PriorRef(p, t)
end

function Pigeons.sample_iid!(target::Comrade.TransformedVLBIPosterior, replica, shared)
    return replica.state = Pigeons.initialization(target, replica.rng, replica.replica_index)
end

# Cube reference: iid draws are just uniforms on [0,1]^n. Every other latent space (flat,
# StdNormal, …) has the pulled-back prior as its reference, so draw from the prior and pull
# it back into the latent space.
function Pigeons.sample_iid!(::PriorRef{P, <:CubeTransport}, replica, shared) where {P}
    return rand!(replica.rng, replica.state)
end

function Pigeons.sample_iid!(target::PriorRef, replica, shared)
    Comrade.PT.latent_pback!(replica.state, target.transform, rand(replica.rng, target.model))
    return replica.state
end


function Pigeons.sample_array(tpost::Comrade.TransformedVLBIPosterior, pt::Pigeons.PT)
    samples = sample_array(pt)
    tbl = mapreduce(hcat, eachslice(samples, dims = (3,), drop = true)) do arr
        s = map(x -> Comrade.transform(tpost, @view(x[begin:(end - 1)])), eachrow(arr))
        return s
    end

    sts = (logdensity = samples[:, end, :] |> vec,)

    return Comrade.PosteriorSamples(tbl, sts; metadata = Dict(:sampler => :Pigeons, :post => tpost))
end


LogDensityProblems.dimension(t::PriorRef) = Comrade.dimension(t.transform)
LogDensityProblems.logdensity(t::PriorRef, x) = t(x)
LogDensityProblems.capabilities(::Type{<:PriorRef}) = LogDensityProblems.LogDensityOrder{0}()

Pigeons.LogDensityProblemsAD.ADgradient(kind::Val, log_potential::Comrade.AbstractVLBIPosterior, replica::Pigeons.Replica) =
    Pigeons.BufferedAD(log_potential, replica.recorders.buffers, Ref(0.0), Ref{Cstring}())


function LogDensityProblems.logdensity_and_gradient(log_potential::Pigeons.BufferedAD{<:Comrade.AbstractVLBIPosterior}, x)
    m = log_potential.enclosed
    b = log_potential.buffer
    ∂ℓ_∂x = fill!(b, zero(eltype(b))) # NB: Enzyme gives erroneous answer if buffer is not zeroed first
    mode = EnzymeCore.WithPrimal(Comrade.admode(m))
    _, y = EnzymeCore.autodiff(
        mode, LogDensityProblems.logdensity, EnzymeCore.Active,
        EnzymeCore.Const(m), EnzymeCore.Duplicated(x, ∂ℓ_∂x)
    )
    return y, ∂ℓ_∂x
end

end
