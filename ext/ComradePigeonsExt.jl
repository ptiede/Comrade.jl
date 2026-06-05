module ComradePigeonsExt

using Comrade

using Pigeons
using AbstractMCMC
using EnzymeCore
using LogDensityProblems
using Random

# Latent-space tags for the wrapped `TransportedDistribution`: `StdFlat` (unconstrained
# ℝⁿ, `asflat`) has `stop === nothing`; the unit hypercube (`ascube`) carries a `StdUniform`.
const FlatTransport = Comrade.TransportedDistribution{<:Any, <:Any, Nothing}
const CubeTransport = Comrade.TransportedDistribution{<:Any, <:Any, <:Comrade.StdUniform}


Pigeons.initialization(tpost::Comrade.TransformedVLBIPosterior, rng::Random.AbstractRNG, ::Int) = prior_sample(rng, tpost)

struct PriorRef{P, T}
    model::P
    transform::T
end

# The Pigeons reference is the prior. `transport_and_logdensity` returns the pulled-back
# *prior* log-density directly: `logpdf(prior, y) + logjac` in the flat space and
# `logpdf(StdUniform, x)` (0 inside the cube, -Inf outside) in the cube space — so a single
# expression covers both references.
function (p::PriorRef)(x)
    return last(Comrade.transport_and_logdensity(p.transform, x))
end

Pigeons.default_explorer(::Comrade.TransformedVLBIPosterior{P, <:CubeTransport}) where {P} =
    SliceSampler()

Pigeons.default_explorer(::Comrade.TransformedVLBIPosterior{P, <:FlatTransport}) where {P} =
    Pigeons.AutoMALA(; default_autodiff_backend = :Enzyme)

function Pigeons.default_reference(tpost::Comrade.TransformedVLBIPosterior)
    t = tpost.transform
    p = tpost.lpost.prior
    return PriorRef(p, t)
end

function Pigeons.sample_iid!(target::Comrade.TransformedVLBIPosterior, replica, shared)
    return replica.state = Pigeons.initialization(target, replica.rng, replica.replica_index)
end

function Pigeons.sample_iid!(target::PriorRef{P, <:FlatTransport}, replica, shared) where {P}
    return replica.state .= Comrade.inverse(target.transform, rand(replica.rng, target.model))
end

function Pigeons.sample_iid!(::PriorRef{P, <:CubeTransport}, replica, shared) where {P}
    return rand!(replica.rng, replica.state)
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
