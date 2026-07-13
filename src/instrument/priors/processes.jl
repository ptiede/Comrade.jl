export AbstractGaussMarkovProcess, OrnsteinUhlenbeck

"""
    AbstractGaussMarkovProcess

Abstract type for continuous-time, stationary Gauss-Markov processes in state-space form.
These are used to build time-correlated instrument priors (see [`GaussMarkovSitePrior`](@ref))
whose log-density factorizes over the time chain, i.e. is `O(n)` in the number of times.

# Interface

A concrete process `p` must implement

  - `stationary_moments(p) -> (μ, P)`: the stationary mean and variance.
  - `transition_moments(p, Δt) -> (Φ, Q)`: the exact discretization over a gap `Δt`
    (in hours), so that `x(t + Δt) | x(t) ~ N(μ + Φ*(x(t) - μ), Q)`.

Hyperparameter handling is generic: any field of the process that is a
`Distributions.Distribution` is treated as a *fitted* hyperparameter with that
distribution as its prior (see [`hyperprior`](@ref Comrade.hyperprior)), and is substituted
with a sampled value via [`materialize`](@ref Comrade.materialize) before
`stationary_moments`/`transition_moments` are called. Fields that are plain numbers are
fixed. The generic `materialize` uses `ConstructionBase.setproperties`, so concrete types
work out of the box and can customize reconstruction via `ConstructionBase.constructorof`.

!!! note
    The interface is written in state-space form so that higher-order processes
    (e.g. Matérn-3/2, CARMA) with vector-valued latent states can be added in the
    future. Vector-state processes will additionally require an `observation(p)`
    method returning the observation vector `H`; this is not implemented yet and
    only scalar-state processes are currently supported.
"""
abstract type AbstractGaussMarkovProcess end

"""
    stationary_moments(p::AbstractGaussMarkovProcess) -> (μ, P)

Return the stationary mean `μ` and variance `P` of the process `p`.
"""
function stationary_moments end

"""
    transition_moments(p::AbstractGaussMarkovProcess, Δt) -> (Φ, Q)

Return the transition coefficient `Φ` and innovation variance `Q` for the exact
discretization of `p` over a time gap `Δt` in hours, so that
`x(t + Δt) | x(t) ~ N(μ + Φ*(x(t) - μ), Q)`.
"""
function transition_moments end

"""
    hyperprior(p::AbstractGaussMarkovProcess) -> NamedTuple

Return the priors of the fitted hyperparameters of `p`, i.e. the fields of `p` that are
`Distributions.Distribution`s. Returns an empty NamedTuple when every field is a number,
in which case the process has no fitted hyperparameters.
"""
function hyperprior(p::AbstractGaussMarkovProcess)
    return filter(Base.Fix2(isa, Dists.Distribution), getproperties(p))
end

"""
    materialize(p::AbstractGaussMarkovProcess, hp::NamedTuple)

Substitute the fitted hyperparameters of `p` with the values in `hp`, returning a fully
numeric process. Fields of `p` that are `Distributions.Distribution`s are replaced by
`hp.<fieldname>`; numeric fields are kept. When `p` has no fitted hyperparameters this
is the identity.
"""
materialize(p::AbstractGaussMarkovProcess, ::NamedTuple{()}) = p
# Subset `hp` to the fitted fields first: chains reading the default hyperparameters are
# handed the whole merged NamedTuple, which may also carry per-site override entries.
materialize(p::AbstractGaussMarkovProcess, hp::NamedTuple) = setproperties(p, hp[keys(hyperprior(p))])

"""
    OrnsteinUhlenbeck(; σ, τ, μ = 0.0)

A stationary Ornstein-Uhlenbeck process, i.e. the continuous-time AR(1) process. This is
the mean-reverting Gaussian process with marginal `x(t) ~ N(μ, σ²)` and exponential
correlation `cov(x(t), x(t+Δ)) = σ² exp(-|Δ|/τ)`, whose exact discretization over a gap
`Δt` is

    x(t + Δt) | x(t) ~ N(μ + φ*(x(t) - μ), σ²*(1 - φ²)),   φ = exp(-Δt/τ).

# Arguments

  - `σ`: the marginal (stationary) standard deviation. Pass a `Number` to fix it or a
    `Distributions.Distribution` to fit it as a hyperparameter with that prior.
  - `τ`: the correlation timescale in **hours** (the times of a VLBI observation are in
    UTC hours). A `Number` fixes it; a `Distributions.Distribution` fits it.
  - `μ`: the stationary mean (default `0.0`). Not fittable currently.

# Example

```julia
# fixed hyperparameters
OrnsteinUhlenbeck(σ = 0.1, τ = 2.0)
# fitted marginal std and correlation time
OrnsteinUhlenbeck(σ = Exponential(0.1), τ = InverseGamma(3.0, 6.0))
```
"""
struct OrnsteinUhlenbeck{Tσ, Tτ, Tμ} <: AbstractGaussMarkovProcess
    σ::Tσ
    τ::Tτ
    μ::Tμ
end

function OrnsteinUhlenbeck(; σ, τ, μ = 0.0)
    σ isa Number && σ ≤ 0 && throw(ArgumentError("OrnsteinUhlenbeck σ must be positive"))
    τ isa Number && τ ≤ 0 && throw(ArgumentError("OrnsteinUhlenbeck τ must be positive"))
    μ isa Dists.Distribution && throw(ArgumentError("OrnsteinUhlenbeck μ is not fittable; pass a number"))
    return OrnsteinUhlenbeck(σ, τ, μ)
end

@inline stationary_moments(p::OrnsteinUhlenbeck) = (p.μ, p.σ^2)

@inline function transition_moments(p::OrnsteinUhlenbeck, Δt)
    Φ = exp(-Δt / p.τ)
    return Φ, p.σ^2 * (1 - Φ^2)
end
