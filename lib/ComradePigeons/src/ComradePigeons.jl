module ComradePigeons

using Reexport
@reexport using Pigeons
using Comrade
using Random
using LogDensityProblemsAD

export PigeonsLogPotential

struct PigeonsLogPotential{M}
    post::M
end

"""
    PigeonsLogPotential(post::M)
    PigeonsLogPotential(post::M, ad_backend::Symbol)

Creates a pigeons logdensity object from a Comrade posterior object that satisfies the Pigeons interface.

# Arguments
  - `post`: A Comrade [`Posterior`](@ref) or [`TransformedPosterior`](@ref) object.
  - `ad_backend`: A symbol denoting the type of autodiff backend. This uses the [LogDensityProblemsAD](https://github.com/tpapp/LogDensityProblemsAD.jl) api.
"""
function PigeonsLogPotential(post::M, ad_backend::Symbol) where {M}
    return PigeonsLogPotential(ADgradient(Val(ad_backend)), post)
end


# This one takes in the log jacobian of the transformation not the prior!
function (m::PigeonsLogPotential)(x)
    return LogDensityProblems.logdensity(m.post, x)
end

function Pigeons.gradient!!(log_potential::PigeonsLogPotential, x::T, buffer::T) where {T}
    buffer .= LogDensityProblemsAD.logdensity_and_gradient(log_potential, x)[2]
    return buffer
end

# Pigeons.@provides target PigeonsLogPotential(model::Comrade.TransformedPosterior) =
#     PigeonsLogPotential(model)

Pigeons.create_state_initializer(target::PigeonsLogPotential, ::Inputs) = target
function Pigeons.initialization(target::PigeonsLogPotential, rng::Pigeons.SplittableRandom, _::Int64)
   return  Comrade.prior_sample(rng, target.post)
end

Pigeons.default_explorer(::PigeonsLogPotential) = Pigeons.SliceSampler()

Pigeons.create_reference_log_potential(target::PigeonsLogPotential, ::Inputs) = PriorPotential(target.post)

function Pigeons.sample_iid!(target::PigeonsLogPotential, replica)
    replica.state = initialization(target, replica.rng, replica.replica_index)
end




struct PriorPotential{M,T}
    prior::M
    transform::T
    function PriorPotential(post::Comrade.TransformedPosterior)
        t = post.transform
        prior = post.lpost.prior
        return new{typeof(prior), typeof(t)}(prior, t)
    end
end

function (m::PriorPotential)(x)
    y, lj = Comrade.transform_and_logjac(m.transform, x)
    return logdensityof(m.prior, y) + lj
end

function Pigeons.sample_iid!(target::PriorPotential, replica)
    replica.state = Comrade.inverse(target.transform, rand(replica.rng, target.prior))
end


function Pigeons.Inputs(tpost::Comrade.TransformedPosterior; kwargs...)
    return Inputs(;target=PigeonsLogPotential(tpost), kwargs...)
end

end
