module ComradePigeonsExt

using Comrade

using Pigeons
using AbstractMCMC
using ADTypes: AutoEnzyme
using EnzymeCore
using LogDensityProblems
using Random

const CubeTransport = Comrade.CubeTransport


Pigeons.initialization(tpost::Comrade.TransformedVLBIPosterior, rng::Random.AbstractRNG, ::Int) = prior_sample(rng, tpost)

struct PriorRef{P, T, M}
    model::P
    transform::T
    admode::M
end

# The Pigeons reference is the prior: `logpdf(::TransportedDistribution, x)` is the
# pulled-back prior log-density in the flat space and the (transport-free) `StdUniform`
# reference density in the cube space.
(p::PriorRef)(x) = Comrade.Dists.logpdf(p.transform, x)

# The cube ([0,1]^n) reference is bounded, so use a slice sampler; every other latent space
# (flat ℝⁿ, StdNormal, …) is unbounded and differentiable, so use gradient-based AutoMALA.
# `P` must be constrained to `VLBIPosterior` (matching the struct's own bound) so this stays
# strictly more specific than the unparametrized method below — an unbounded `where {P}` is
# instead *ambiguous* with it (more specific in the transport, less specific in `P`).
Pigeons.default_explorer(::Comrade.TransformedVLBIPosterior{P, <:CubeTransport}) where {P <: Comrade.VLBIPosterior} =
    SliceSampler()

Pigeons.default_explorer(::Comrade.TransformedVLBIPosterior) =
    Pigeons.AutoMALA(; default_autodiff_backend = :Enzyme)

function Pigeons.default_reference(tpost::Comrade.TransformedVLBIPosterior)
    t = tpost.transform
    p = tpost.lpost.prior
    # Carry the posterior's AD mode so the reference (prior) gradient uses the same
    # runtime-activity Enzyme configuration as the target on the tempered chains.
    return PriorRef(p, t, Comrade.admode(tpost))
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

# --- gradient interface for AutoMALA / gradient-based explorers -------------------------------
# Pigeons translates the AutoMALA backend Symbol (`:Enzyme`) into an `ADTypes.AutoEnzyme` before
# it builds a gradient wrapper for each log-potential (the interpolated potential wraps a target
# and a reference gradient). Its default Enzyme path (`PigeonsEnzymeExt`) hardcodes
# `ReverseWithPrimal`, which lacks the runtime activity Comrade posteriors need — so a raw
# `AutoEnzyme` there errors. Dispatch on the Comrade types (more specific than PigeonsEnzymeExt's
# generic method) to route both the target posterior and the prior reference through Comrade's
# configured Enzyme mode instead. (`AutoEnzyme` comes from `ADTypes`, a co-trigger of this
# extension — `Pigeons.ADTypes` is not reliably reachable across Julia/Pigeons versions.)
Pigeons.LogDensityProblemsAD.ADgradient(::AutoEnzyme, lp::Comrade.AbstractVLBIPosterior, buffers::Pigeons.Augmentation) =
    Pigeons.BufferedAD(lp, buffers)
Pigeons.LogDensityProblemsAD.ADgradient(::AutoEnzyme, lp::PriorRef, buffers::Pigeons.Augmentation) =
    Pigeons.BufferedAD(lp, buffers)

# Reverse-mode gradient of `logdensity(ℓ, x)` into `buffer`, using the given Enzyme `mode` (which
# carries Comrade's runtime-activity setting).
function _enzyme_logdensity_and_gradient(mode, ℓ, x, buffer)
    ∂ℓ_∂x = fill!(buffer, zero(eltype(buffer))) # NB: Enzyme gives an erroneous answer if the buffer is not zeroed first
    _, y = EnzymeCore.autodiff(
        EnzymeCore.WithPrimal(mode), LogDensityProblems.logdensity, EnzymeCore.Active,
        EnzymeCore.Const(ℓ), EnzymeCore.Duplicated(x, ∂ℓ_∂x)
    )
    return y, ∂ℓ_∂x
end

LogDensityProblems.logdensity_and_gradient(b::Pigeons.BufferedAD{<:Comrade.AbstractVLBIPosterior}, x) =
    _enzyme_logdensity_and_gradient(Comrade.admode(b.enclosed), b.enclosed, x, b.buffer)
LogDensityProblems.logdensity_and_gradient(b::Pigeons.BufferedAD{<:PriorRef}, x) =
    _enzyme_logdensity_and_gradient(b.enclosed.admode, b.enclosed, x, b.buffer)

end
