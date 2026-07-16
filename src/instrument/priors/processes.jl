export AbstractGaussMarkovProcess, OrnsteinUhlenbeck, Wiener

"""
    AbstractGaussMarkovProcess

Abstract type for continuous-time, time-homogeneous Gauss-Markov processes in state-space
form. These are the building blocks of time-correlated instrument priors (see
[`GaussMarkovSitePrior`](@ref) and [`GaussMarkovChain`](@ref Comrade.GaussMarkovChain)):
the chain log-density factorizes over the time chain using the process's exact
discretization, so the cost is `O(n)` in the number of times and irregular time sampling
is handled exactly.

# Interface

A concrete process `p` must implement

  - `transition_moments(p, Δt) -> (Φ, Q)`: the exact discretization over a time gap `Δt`
    (in hours), so that `x(t + Δt) | x(t) ~ N(μ + Φ*(x(t) - μ), Q)` with
    `μ = process_mean(p)`.
  - `initial_moments(p) -> (μ₀, P₀)`: the marginal law of the chain's first point.
    Stationary processes get this for free by defining
    [`stationary_moments`](@ref Comrade.stationary_moments), to which `initial_moments`
    falls back.

and optionally

  - `process_mean(p)`: the mean level `μ` the transitions revert to (default `0.0`).
  - `has_proper_initial(p)`: return `false` when the process has no proper initial law
    (e.g. [`Wiener`](@ref)). Such a process requires the chain's first point to be
    exactly conditioned (`anchored = true` or a reference fixing) and `initial_moments`
    is then never called.

The process does *not* need to be stationary: only the Gauss-Markov property and
time-homogeneity are used.

Hyperparameter handling is generic: any field of the process that is a
`Distributions.Distribution` is treated as a *fitted* hyperparameter with that
distribution as its prior (see [`hyperprior`](@ref Comrade.hyperprior)), and is substituted
with a sampled value via [`materialize`](@ref Comrade.materialize) before
the moment functions are called. Fields that are plain numbers are fixed. The generic
`materialize` uses `ConstructionBase.setproperties`, so concrete types work out of the
box and can customize reconstruction via `ConstructionBase.constructorof`.

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

Return the stationary mean `μ` and variance `P` of a *stationary* process `p`.
Non-stationary processes (e.g. [`Wiener`](@ref)) do not define this; the chain machinery
only ever uses [`initial_moments`](@ref Comrade.initial_moments), which falls back to
this function.
"""
function stationary_moments end

"""
    initial_moments(p::AbstractGaussMarkovProcess) -> (μ₀, P₀)

Return the mean and variance of the marginal law of the chain's first point. Defaults to
[`stationary_moments`](@ref Comrade.stationary_moments). Only called when
[`has_proper_initial`](@ref Comrade.has_proper_initial) is `true`.
"""
initial_moments(p::AbstractGaussMarkovProcess) = stationary_moments(p)

"""
    has_proper_initial(p::AbstractGaussMarkovProcess) -> Bool

Whether the process has a proper (finite-variance) initial law. When `false` (e.g.
[`Wiener`](@ref)) the chain's first point must be exactly conditioned — via
`anchored = true` or a reference fixing that covers it — and
[`initial_moments`](@ref Comrade.initial_moments) is never called. Defaults to `true`.
"""
has_proper_initial(::AbstractGaussMarkovProcess) = true

"""
    process_mean(p::AbstractGaussMarkovProcess)

The mean level `μ` the process transitions revert to, i.e. the `μ` in
`x(t + Δt) | x(t) ~ N(μ + Φ*(x(t) - μ), Q)`. Defaults to `0.0`; irrelevant for
processes with `Φ = 1` (e.g. [`Wiener`](@ref)).
"""
process_mean(::AbstractGaussMarkovProcess) = 0.0

"""
    transition_moments(p::AbstractGaussMarkovProcess, Δt) -> (Φ, Q)

Return the transition coefficient `Φ` and innovation variance `Q` for the exact
discretization of `p` over a time gap `Δt` in hours, so that
`x(t + Δt) | x(t) ~ N(μ + Φ*(x(t) - μ), Q)` with `μ = process_mean(p)`.
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
# Subset `hp` to the fitted fields so extra entries are ignored rather than an error.
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
@inline process_mean(p::OrnsteinUhlenbeck) = p.μ

@inline function transition_moments(p::OrnsteinUhlenbeck, Δt)
    # Keep Φ strictly below 1 so that Q and the conditional (bridge) variances downstream
    # stay positive even when Δt/τ underflows, e.g. during a fitted-τ warmup excursion.
    # Otherwise the chain log-density is NaN instead of merely very negative, and the
    # whitening transform divides by zero.
    Φ = min(exp(-Δt / p.τ), prevfloat(1.0))
    return Φ, p.σ^2 * (1 - Φ^2)
end

"""
    Wiener(; σ)

A Wiener process (Brownian motion), i.e. the continuous-time random walk with
`x(t + Δt) | x(t) ~ N(x(t), σ² Δt)` and structure function
`E[(x(t+Δ) - x(t))²] = σ² Δ` — the classical model for atmospheric phase drift.

The process is non-stationary and has no proper initial law, so a chain using it must
have its first point exactly conditioned: use
`GaussMarkovSitePrior(seg, Wiener(σ = ...); anchored = true)` (or a reference-fixing
scheme that covers each chain's first point). The prior then describes the evolution of
the quantity relative to its starting value.

# Arguments

  - `σ`: the diffusion scale per √hour (the times of a VLBI observation are in UTC
    hours). Pass a `Number` to fix it or a `Distributions.Distribution` to fit it as a
    hyperparameter with that prior.

# Example

```julia
dgp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), Wiener(σ = 1.0); anchored = true); refant = SEFDReference(0.0))
```
"""
struct Wiener{Tσ} <: AbstractGaussMarkovProcess
    σ::Tσ
end

function Wiener(; σ)
    σ isa Number && σ ≤ 0 && throw(ArgumentError("Wiener σ must be positive"))
    return Wiener(σ)
end

has_proper_initial(::Wiener) = false

@inline transition_moments(p::Wiener, Δt) = (one(Δt), p.σ^2 * Δt)
