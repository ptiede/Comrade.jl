export AbstractGaussMarkovProcess, OrnsteinUhlenbeck, BrownianMotion, WrappedBrownian,
    WrappedOrnsteinUhlenbeck, VonMisesProcess
export AbstractInitialPrior, StationaryInit, FixedInit, GaussianInit, UniformInit

"""
    AbstractGaussMarkovProcess

Abstract type for continuous-time Gauss-Markov processes in state-space form, used to
build time-correlated instrument priors (see [`GaussMarkovSitePrior`](@ref)) whose
log-density factorizes over the time chain and so costs `O(n)` in the number of times.

# Interface

A concrete process `p` must implement

  - `isstationary(p) -> Bool`: whether the process has a stationary distribution.
  - `stationary_moments(p) -> (μ, P)`: the stationary mean and variance (stationary
    processes only).
  - `transition_moments(p, Δt) -> (Φ, Q)`: the exact discretization over a gap `Δt`
    (in hours), so that `x(t + Δt) | x(t) ~ N(μ + Φ*(x(t) - μ), Q)`.

The distribution of each chain's first time stamp is set by the site prior's
[`AbstractInitialPrior`](@ref), not by the process.

Any field of the process that is a `Distributions.Distribution` is a fitted hyperparameter
with that distribution as its prior (see [`hyperprior`](@ref Comrade.hyperprior)); it is
replaced by a sampled value via [`materialize`](@ref Comrade.materialize) before
`stationary_moments`/`transition_moments` are called. Numeric fields are fixed. The generic
`materialize` uses `ConstructionBase.setproperties`; a concrete type can customize
reconstruction via `ConstructionBase.constructorof`.

!!! note
    Only scalar-state processes are supported.
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
    isstationary(p::AbstractGaussMarkovProcess) -> Bool

Whether `p` has a stationary distribution, i.e. implements
[`stationary_moments`](@ref Comrade.stationary_moments). [`StationaryInit`](@ref), the
default initial prior of [`GaussMarkovSitePrior`](@ref), requires it. Every process must
implement this; there is no fallback.
"""
function isstationary end

# Whether the process lives on the circle. A wrapped process's moments are read as wrapped
# normal (see `_marginal_logpdf`), and it supports only `asflat`.
is_wrapped(::AbstractGaussMarkovProcess) = false

# Whether a reference-fixed point is a gauge choice rather than a value of the process
# (see "Reference stations" in the `GaussMarkovSitePrior` docstring). For `WrappedBrownian`
# the free-stamp chain is the exact marginal of the full one; for a mean-reverting wrapped
# process it is not, since the shortest-arc drift does not compose.
ref_is_gauge(p::AbstractGaussMarkovProcess) = is_wrapped(p)

# Whether the non-centered lift resolves with `affine_scan!`, which needs a transition
# mean affine in the previous point: every real-line process, and `WrappedBrownian`
# (Φ ≡ 1). `VonMisesProcess` overrides this to false.
has_affine_lift(p::AbstractGaussMarkovProcess) = !is_wrapped(p) || !_periodic_mean(p)

# Whether the process supports the non-centered lift, which needs a transition mean that
# carries the previous point's sheet: affine in it, or `2π`-equivariant
# (`m(xp + 2π) = m(xp) + 2π`) like `VonMisesProcess`. The shortest-arc drift of
# `WrappedOrnsteinUhlenbeck` (`_periodic_mean`) does not.
has_noncentered_lift(p::AbstractGaussMarkovProcess) = has_affine_lift(p)

# Mean entering the transition law `x(t+Δ) | x(t) ~ N(μ + Φ(x(t)−μ), Q)`. Stationary
# processes revert to their stationary mean. Nonstationary processes (Φ ≡ 1) must
# override with any finite value; it cancels out of every Φ = 1 formula.
_process_mean(p::AbstractGaussMarkovProcess) = first(stationary_moments(p))

# Transition mean `μ + Φ(xp − μ)` given the previous point; wrapped mean-reverting
# processes override it.
@inline _cond_mean(::AbstractGaussMarkovProcess, xp, Φ, μ) = μ + Φ * (xp - μ)

# Whether `_cond_mean` is `2π`-periodic in the previous point (true for the shortest-arc
# drift). Decides whether the centered lift's inverse is pointwise or a recurrence.
_periodic_mean(::AbstractGaussMarkovProcess) = false

# Variance floor for the chain log-densities. A zero variance (a fitted hyperparameter at
# the edge of its support, e.g. a τ quantile of `Inf`, or underflow) makes `_normal_logpdf`
# and `_wn_logpdf` return NaN; floored, the density is finite and very negative, so the
# sampler rejects the point. `eps²` is below any variance a positive hyperparameter
# produces (the OU Φ clamp gives at least `σ²·eps`).
@inline _floor_var(Q) = max(Q, eps(one(Q))^2)

# Log-density of a deviation `δ` with variance `P`: normal on the real line, wrapped normal
# on the circle. The `Val` dispatch constant-folds, keeping traced loop bodies single-path.
@inline _marginal_logpdf(p::AbstractGaussMarkovProcess, δ, P) =
    _marginal_logpdf(Val(is_wrapped(p)), δ, P)
@inline _marginal_logpdf(::Val{false}, δ, P) = _normal_logpdf(δ, P)
@inline _marginal_logpdf(::Val{true}, δ, P) = _wn_logpdf(δ, P)

# Transition log-density over a gap Δt, from the process's own transition moments.
@inline function _transition_logpdf(p::AbstractGaussMarkovProcess, xi, xp, Δt, μ)
    Φ, Q = transition_moments(p, Δt)
    return _marginal_logpdf(p, xi - _cond_mean(p, xp, Φ, μ), Q)
end

"""
    hyperprior(p::AbstractGaussMarkovProcess) -> NamedTuple

Return the priors of the fitted hyperparameters of `p`, i.e. the fields of `p` that are
`Distributions.Distribution`s. Returns an empty NamedTuple when every field is a number,
in which case the process has no fitted hyperparameters.
"""
function hyperprior(p::AbstractGaussMarkovProcess)
    return filter(Base.Fix2(isa, Dists.Distribution), getproperties(p))
end

# Float type of a materialized sample of `p`: numeric fields contribute their type, fitted
# fields their `eltype`.
_paramtype(p::AbstractGaussMarkovProcess) =
    float(mapreduce(_ptype, promote_type, values(getproperties(p)); init = Union{}))
_ptype(x::Number) = typeof(x)
_ptype(d::Dists.Distribution) = eltype(d)

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

# A fixed hyperparameter must be positive and finite (NaN fails the comparison too): zero
# or Inf puts a zero in a transition variance. Fitted fields are guarded by `_floor_var`.
function _check_hyper(process::Symbol, name::Symbol, x)
    x isa Number && !(0 < x < Inf) &&
        throw(ArgumentError("$process $name must be positive and finite"))
    return nothing
end

# The mean is not fittable and enters the transition mean linearly, so it need only be
# a finite number.
function _check_mean(process::Symbol, μ)
    μ isa Dists.Distribution &&
        throw(ArgumentError("$process μ is not fittable; pass a number"))
    isfinite(μ) || throw(ArgumentError("$process μ must be finite"))
    return nothing
end

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
  - `μ`: the stationary mean (default `0.0`). Not fittable.

!!! note
    Fitted hyperparameters are evaluated pointwise, so on GPU/`Reactant` backends their
    priors must accept a scalar traced argument. Use the `VLBI*` distributions
    (e.g. `VLBIExponential`, `VLBIInverseGamma`) rather than the plain `Distributions.jl`
    types for a model you intend to run through `Reactant`.

# Example

```julia
# fixed hyperparameters
OrnsteinUhlenbeck(σ = 0.1, τ = 2.0)
# fitted marginal std and correlation time (VLBI* variants work on CPU and Reactant)
OrnsteinUhlenbeck(σ = VLBIExponential(0.1), τ = VLBIInverseGamma(3.0, 6.0))
```
"""
struct OrnsteinUhlenbeck{Tσ, Tτ, Tμ} <: AbstractGaussMarkovProcess
    σ::Tσ
    τ::Tτ
    μ::Tμ
end

function OrnsteinUhlenbeck(; σ, τ, μ = 0.0)
    _check_hyper(:OrnsteinUhlenbeck, :σ, σ)
    _check_hyper(:OrnsteinUhlenbeck, :τ, τ)
    _check_mean(:OrnsteinUhlenbeck, μ)
    return OrnsteinUhlenbeck(σ, τ, μ)
end

isstationary(::OrnsteinUhlenbeck) = true

@inline stationary_moments(p::OrnsteinUhlenbeck) = (p.μ, p.σ^2)

# Exact OU discretization, shared by the OU-family processes. Φ must stay below 1 so Q and
# the bridge variances stay positive when Δt/τ underflows (e.g. a very large fitted τ);
# at Φ = 1 the log-density is NaN and the whitening divides by zero.
@inline function _ou_transition_moments(σ, τ, Δt)
    Φ = min(exp(-Δt / τ), prevfloat(1.0))
    return Φ, σ^2 * (1 - Φ^2)
end

@inline transition_moments(p::OrnsteinUhlenbeck, Δt) = _ou_transition_moments(p.σ, p.τ, Δt)

"""
    BrownianMotion(; D)

A driftless Brownian motion (Wiener process) with diffusion coefficient `D`: increments
over a gap `Δt` (in hours) are independent `N(0, D*Δt)`, so `√D` is the typical drift per
`√hour`. This is the non-reverting limit of [`OrnsteinUhlenbeck`](@ref) (`τ → ∞` at fixed
`D = 2σ²/τ`), for quantities that wander without a preferred level.

The process is nonstationary (the marginal variance grows linearly in time), so a chain
needs an explicit initial prior, [`GaussianInit`](@ref) or [`FixedInit`](@ref);
`StationaryInit` throws.

# Arguments

  - `D`: the diffusion coefficient, in units of the parameter squared per hour. Pass a
    `Number` to fix it or a `Distributions.Distribution` to fit it as a hyperparameter
    with that prior.
"""
struct BrownianMotion{TD} <: AbstractGaussMarkovProcess
    D::TD
end

function BrownianMotion(; D)
    _check_hyper(:BrownianMotion, :D, D)
    return BrownianMotion(D)
end

isstationary(::BrownianMotion) = false
stationary_moments(::BrownianMotion) =
    throw(ArgumentError("BrownianMotion has no stationary distribution"))
_process_mean(p::BrownianMotion) = zero(_paramtype(p))

@inline transition_moments(p::BrownianMotion, Δt) = (one(Δt), p.D * Δt)

"""
    WrappedBrownian(; τ)

Brownian motion on the circle: angular increments over a gap `Δt` (hours) are wrapped
normal `WN(0, 2Δt/τ)`, the wrap of [`BrownianMotion`](@ref) onto `(−π, π]`. Use it for
gain phases: the prior is `2π`-periodic, so `2π`-shifted posterior modes are equivalent
(the likelihood sees only `exp(iθ)`), and the coherence-time posterior is not biased by
phase wraps.

In the visibility domain the process decorrelates exponentially,
`E[e^{i(θ(t+Δ) − θ(t))}] = e^{−Δ/τ}`, so `τ` is the phase coherence time in hours (the
diffusion coefficient is `D = 2/τ` rad² per hour). The circular stationary law is
uniform, so with [`UniformInit`](@ref) the chain also absorbs the per-track phase offset
and no separate circular offset term is needed.

The transition mean is the previous phase, so this process supports
`param = NonCentered()` (the default): each free phase's coordinate is its increment
divided by the step's standard deviation, and the phase is the wrapped running sum.
`param = Centered()` uses the raw angles. Both use sheet weights to keep the flat target
proper without changing the circular model; [`AngleEmbedded`](@ref) avoids them at twice
the phase dimension. Only `asflat` is supported. See [`GaussMarkovSitePrior`](@ref).

A reference-fixed time stamp is a gauge choice, so it is dropped from the chain and the
free points step over the combined gap (wrapped normal transitions compose: their
variances add). A reference station's phases at other times are offsets from its last
free phase, not pulled toward the reference value.

# Arguments

  - `τ`: the phase coherence time in **hours** (the times of a VLBI observation are in
    UTC hours): the gap over which the visibility-domain coherence
    `E[e^{iΔθ}]` falls by `1/e`. Pass a `Number` to fix it or a
    `Distributions.Distribution` to fit it as a hyperparameter with that prior.

# Example

```julia
gp ~ ArrayPrior(
    GaussMarkovSitePrior(IntegSeg(), WrappedBrownian(τ = VLBIInverseGamma(3.0, 6.0)); init = UniformInit());
    refant = SEFDReference(0.0)
)
```
"""
struct WrappedBrownian{Tτ} <: AbstractGaussMarkovProcess
    τ::Tτ
end

function WrappedBrownian(; τ)
    _check_hyper(:WrappedBrownian, :τ, τ)
    return WrappedBrownian(τ)
end

isstationary(::WrappedBrownian) = false
is_wrapped(::WrappedBrownian) = true
stationary_moments(::WrappedBrownian) = throw(
    ArgumentError(
        "WrappedBrownian has no Gaussian stationary moments; its circular stationary " *
            "law is uniform (UniformInit)"
    )
)
_process_mean(p::WrappedBrownian) = zero(_paramtype(p))
# Increment variance D*Δt with D = 2/τ, so the coherence exp(-Q/2) e-folds at Δt = τ.
# Floored because Δt/τ can underflow, and the conditional variances and `√Q` must stay
# positive.
@inline _increment_var(p::WrappedBrownian, Δt) = _floor_var(2 * Δt / p.τ)
@inline transition_moments(p::WrappedBrownian, Δt) = (one(Δt), _increment_var(p, Δt))

"""
    WrappedOrnsteinUhlenbeck(; σ, τ, μ = 0.0)

The mean-reverting process on the circle: the wrap of [`OrnsteinUhlenbeck`](@ref) onto
`(−π, π]`, with the drift taken along the **shortest arc** so that the law is
`2π`-periodic in both of its arguments and is a circular Markov kernel,

    θ(t + Δt) | θ(t) ~ WN(μ + φ*wrap(θ(t) - μ), σ²*(1 - φ²)),   φ = exp(-Δt/τ),

where `wrap` maps to `(−π, π]`. Its circular marginal is the wrapped normal `WN(μ, σ²)`,
so, unlike [`WrappedBrownian`](@ref), the phase stays near a level, and the
visibility-domain coherence

    E[e^{i(θ(t+Δ) − θ(t))}] = exp(−σ²(1 − e^{−Δ/τ}))

levels off at `exp(−σ²)` rather than decaying to zero. Taking `σ, τ → ∞` at fixed
`2σ²/τ = 2/τ_coh` recovers `WrappedBrownian(τ = τ_coh)`.

Use it for a gain ratio (R/L) phase, an instrumental offset that is stable up to slow
drift: one chain replaces a [`WrappedBrownian`](@ref) chain with [`FixedInit`](@ref) plus a
separate circular offset, and `σ` sets how far the ratio phase may stray. The shortest-arc
drift drops the previous phase's sheet, so this process supports only `param = Centered()`
and `asflat`; [`VonMisesProcess`](@ref) is a non-centered alternative. As for every
wrapped process, reference-fixed time stamps are dropped from the chain. See
[`GaussMarkovSitePrior`](@ref).

# Arguments

  - `σ`: the circular marginal spread in **radians**, i.e. the `WN(μ, σ²)` scale. Pass a
    `Number` to fix it or a `Distributions.Distribution` to fit it as a hyperparameter.
  - `τ`: the reversion timescale in **hours** (the times of a VLBI observation are in UTC
    hours), the same units as the `τ` of [`OrnsteinUhlenbeck`](@ref) and
    [`WrappedBrownian`](@ref). A `Number` fixes it; a `Distributions.Distribution` fits it.
  - `μ`: the circular mean (default `0.0`). Not fittable.

!!! warning "Valid for `σ` well below `π`"
    Three properties hold only up to the `WN(μ, σ²)` mass beyond `μ ± π`: stationarity of
    the marginal, exact composition of the reference-antenna conditioning of
    [`GaussMarkovSitePrior`](@ref) over skipped time stamps, and continuity of the kernel
    at the antipode of `μ`. That mass is `1.7e-3` at `σ = 1` rad and `e^{-123}` at
    `σ = 0.2` rad. Do not use this process as an uninformative circular prior (`σ ≳ 2`);
    use [`WrappedBrownian`](@ref) with [`UniformInit`](@ref) instead.

# Example

```julia
## an R/L ratio phase: pinned near zero, drifting on a 12 hour timescale
gprat ~ ArrayPrior(
    GaussMarkovSitePrior(
        ScanSeg(), WrappedOrnsteinUhlenbeck(σ = VLBIExponential(0.2), τ = VLBIInverseGamma(1.0, 12.0))
    );
    refant = SingleReference(:AA, 0.0)
)
```
"""
struct WrappedOrnsteinUhlenbeck{Tσ, Tτ, Tμ} <: AbstractGaussMarkovProcess
    σ::Tσ
    τ::Tτ
    μ::Tμ
end

function WrappedOrnsteinUhlenbeck(; σ, τ, μ = 0.0)
    _check_hyper(:WrappedOrnsteinUhlenbeck, :σ, σ)
    _check_hyper(:WrappedOrnsteinUhlenbeck, :τ, τ)
    _check_mean(:WrappedOrnsteinUhlenbeck, μ)
    return WrappedOrnsteinUhlenbeck(σ, τ, μ)
end

isstationary(::WrappedOrnsteinUhlenbeck) = true
is_wrapped(::WrappedOrnsteinUhlenbeck) = true
# Read as the *wrapped* normal WN(μ, σ²) by `_marginal_logpdf`, since `is_wrapped` holds.
@inline stationary_moments(p::WrappedOrnsteinUhlenbeck) = (p.μ, p.σ^2)

# The same exact discretization as `OrnsteinUhlenbeck`, read as wrapped normal.
@inline transition_moments(p::WrappedOrnsteinUhlenbeck, Δt) = _ou_transition_moments(p.σ, p.τ, Δt)

# The shortest-arc drift makes the kernel 2π-periodic in the previous angle.
@inline _cond_mean(p::WrappedOrnsteinUhlenbeck, xp, Φ, μ) = μ + Φ * _wrap_angle(xp - μ)
_periodic_mean(::WrappedOrnsteinUhlenbeck) = true

"""
    VonMisesProcess(; σ, τ, μ = 0.0)

The mean-reverting circular process with a smooth drift. It has the same `σ`, `τ`, `μ`,
and transition variance `σ²(1 − φ²)`, `φ = exp(−Δt/τ)`, as
[`WrappedOrnsteinUhlenbeck`](@ref), with the shortest-arc drift replaced by the sine drift

    θ(t + Δt) | θ(t) ~ WN(θ(t) − (1 − φ)·sin(θ(t) − μ), σ²*(1 − φ²)),

a discretization of the von Mises (circular Langevin) diffusion, whose stationary law is
`vM(μ, 1/σ²) ≈ WN(μ, σ²)` for small `σ`. For deviations well inside `(−π, π)` the two
drifts agree (`sin ε ≈ ε`), so this process can replace `WrappedOrnsteinUhlenbeck` in the
`σ ≪ π` regime, with the same warning: do not use it as an uninformative circular prior;
use [`WrappedBrownian`](@ref) there.

The sine mean is smooth and `2π`-equivariant in the previous point, whereas the
shortest-arc drift jumps by `2πφ` at the antipode of `μ`, so this process supports
`param = NonCentered()` (the default), whose coordinates decouple the chain values from
fitted hyperparameters and from each other. The drift also vanishes at the antipode, an
unstable equilibrium where the stationary density is smallest; for `σ ≪ π` the chain
almost never reaches it.

The drift is nonlinear, so a non-centered chain is resolved as a sequential recurrence
rather than the batched affine scan; under `Reactant` this is a serialized loop. It suits
scan-level chains, e.g. a gain ratio phase about a separate circular offset with
`init = FixedInit(0.0)`, or carrying its own level with the default `StationaryInit`. For
an integration-level absolute phase use [`WrappedBrownian`](@ref).

# Arguments

  - `σ`: the circular marginal spread in **radians**. Pass a `Number` to fix it or a
    `Distributions.Distribution` to fit it as a hyperparameter.
  - `τ`: the reversion timescale in **hours** (the times of a VLBI observation are in UTC
    hours). A `Number` fixes it; a `Distributions.Distribution` fits it.
  - `μ`: the circular mean (default `0.0`). Not fittable.

# Example

```julia
## an R/L ratio phase: pinned near zero, drifting on a 12 hour timescale, whitened
gprat ~ ArrayPrior(
    GaussMarkovSitePrior(
        ScanSeg(), VonMisesProcess(σ = VLBIExponential(0.2), τ = VLBIInverseGamma(1.0, 12.0))
    );
    refant = SingleReference(:AA, 0.0)
)
```
"""
struct VonMisesProcess{Tσ, Tτ, Tμ} <: AbstractGaussMarkovProcess
    σ::Tσ
    τ::Tτ
    μ::Tμ
end

function VonMisesProcess(; σ, τ, μ = 0.0)
    _check_hyper(:VonMisesProcess, :σ, σ)
    _check_hyper(:VonMisesProcess, :τ, τ)
    _check_mean(:VonMisesProcess, μ)
    return VonMisesProcess(σ, τ, μ)
end

isstationary(::VonMisesProcess) = true
is_wrapped(::VonMisesProcess) = true
# Read as WN(μ, σ²), the small-σ form of the vM(μ, 1/σ²) stationary law.
@inline stationary_moments(p::VonMisesProcess) = (p.μ, p.σ^2)
@inline transition_moments(p::VonMisesProcess, Δt) = _ou_transition_moments(p.σ, p.τ, Δt)

# The sine drift: smooth and 2π-equivariant in the previous point (m(xp + 2π) =
# m(xp) + 2π), matching the linear OU drift μ + Φ(xp − μ) to first order about μ.
# Monotone in xp (m′ = 1 − (1 − Φ)cos ≥ Φ > 0), so the mean map is a circle
# diffeomorphism and the non-centered recurrence never folds.
@inline _cond_mean(p::VonMisesProcess, xp, Φ, μ) = xp - (1 - Φ) * sin(xp - μ)
# Equivariant but not periodic, so the non-centered lift exists; not affine, so it
# resolves sequentially.
has_affine_lift(::VonMisesProcess) = false
has_noncentered_lift(::VonMisesProcess) = true

# Wrap an angle to (−π, π] up to floating point; `round` keeps the map traceable and its
# a.e. zero derivative makes the wrap transparent to AD.
@inline _wrap_angle(θ) = θ - oftype(float(θ), 2π) * round(θ / oftype(float(θ), 2π))

# Wrapped-normal log-density of an angular deviation Δ with variance Q: the image sum for
# small Q, the Fourier sum for large Q, each exact to machine precision at the crossover
# Q = 4 with three terms. Both are evaluated and selected with `ifelse` to keep one trace
# under Reactant; the Fourier sum can go negative for small Q, where it is not selected,
# so it is floored before the log. The image sum factors out its k = 0 term (the largest,
# since |Δw| ≤ π), which avoids `log(0)` once abs2(Δw)/(2Q) ≳ 745.
@inline function _wn_logpdf(Δ, Q)
    T2π = oftype(float(Δ), 2π)
    Qf = _floor_var(Q)
    Δw = _wrap_angle(Δ)
    sratio = sum(k -> exp(-(abs2(Δw + T2π * k) - abs2(Δw)) / (2Qf)), -3:3)
    limg = -abs2(Δw) / (2Qf) + log(sratio) - log(T2π * Qf) / 2
    sfouri = ntuple(Val(3)) do j
        Base.@_inline_meta
        exp(-j^2 * Qf / 2) * cos(j * Δw)
    end
    sfour = 1 + 2 * sum(sfouri)
    lfour = log(max(sfour, eps(one(sfour)))) - log(T2π)
    return ifelse(Qf < 4, limg, lfour)
end

# ----- sheet weights: proper flat lift of a wrapped chain ---------------------
#
# A wrapped chain's density is `2π`-periodic in every free angle, so reading a real
# coordinate directly as the angle gives an improper flat target with identical modes
# spaced `2π` apart. The flat density of a free coordinate `φ` (angle `φ mod 2π`) is
# therefore multiplied by
#
#     w(φ) = N(φ − ref; 0, Q) / WN(φ − ref; Q)
#
# for a reference `ref` that does not depend on the sheet of `φ`. Since
# `Σₘ N(δ + 2πm; 0, Q) = WN(δ; Q)`, `Σₘ w(φ + 2πm) = 1`: the pushforward under
# `φ ↦ φ mod 2π` is the chain density, and `w` only concentrates the mass on the sheet
# nearest `ref`. With `ref` the previous chain point and `Q` its transition variance, the
# weight cancels that point's wrapped-normal transition, leaving the unwrapped Gaussian
# random walk `Πₖ N(φₖ − φₖ₋₁; 0, Qₖ)` times the periodic likelihood.

# Log-density of `N(Δ; 0, Q)`. `Q` is floored as in `_wn_logpdf`, so the sheet weights
# still sum to one at the floor.
@inline function _normal_logpdf(Δ, Q)
    Qf = _floor_var(Q)
    return -(abs2(Δ) / Qf + log(oftype(float(Δ), 2π) * Qf)) / 2
end

# A chain's first free point has no earlier value to anchor on; its weight
# (`_init_sheet_logweight`) is defined with the initial distributions below.

# Log sheet weight `log N(δ; 0, Q) − log WN(δ; Q)` of a wrapped chain's free point `xi`,
# anchored on the previous value `xp` a gap `Δt` earlier, with `δ` the deviation from
# `_cond_mean`. Both densities use the same `_cond_mean`, so the weights sum to one over
# the sheets for any process.
@inline function _sheet_logweight(p::AbstractGaussMarkovProcess, xi, xp, Δt, μ)
    Φ, Q = transition_moments(p, Δt)
    return _sheet_logweight_delta(xi - _cond_mean(p, xp, Φ, μ), Q)
end

# The sheet weight on the deviation itself, for the non-centered lift.
@inline _sheet_logweight_delta(δ, Q) = _normal_logpdf(δ, Q) - _wn_logpdf(δ, Q)


# ----- initial distributions p(x(t₁)) ------------------------------------------

"""
    AbstractInitialPrior

Abstract type for the initial distribution `p(x(t₁))` of a [`GaussMarkovSitePrior`](@ref)
chain, i.e. how each (site, frequency) chain's first time stamp is treated:

  - [`StationaryInit`](@ref): the process's stationary marginal (the default; requires a
    stationary process).
  - [`GaussianInit`](@ref): an explicit normal marginal `N(μ0, σ0²)`.
  - [`FixedInit`](@ref): the first point is held at an exact value by conditioning.
  - [`UniformInit`](@ref): uniform on the circle, for nonstationary wrapped processes.

On the real line, initial priors must be Gaussian (a delta is the zero-variance limit):
the exact chain density, the reference-antenna conditioning, and the whitened transports
use closed-form Gaussian conditioning.

# Interface

  - `initial_moments(init, p) -> (m₀, P₀)`: mean and variance of the first time stamp.
"""
abstract type AbstractInitialPrior end

"""
    StationaryInit()

Start each chain in the process's stationary marginal from
[`stationary_moments`](@ref Comrade.stationary_moments): `N(μ, P)` for a real-line process
and the wrapped normal `WN(μ, P)` for a circular one (e.g.
[`WrappedOrnsteinUhlenbeck`](@ref)). Every time stamp then has the same marginal. The
default of [`GaussMarkovSitePrior`](@ref). Requires `isstationary(process)`.
"""
struct StationaryInit <: AbstractInitialPrior end

"""
    FixedInit(value)

Hold each chain at `value` at its first time stamp by exact conditioning, so the prior
describes the quantity's evolution from a known start. Composes with reference-fixed
(`refant`) indices; if the reference scheme fixes a chain's first point to a different
value, posterior construction throws.

`FixedInit(0.0)` on a phase-fluctuation chain on top of a separate circular offset term
removes the level redundancy between the two, which the likelihood's `2π` periodicity
would otherwise turn into posterior modes at `2π`-shifted levels.
"""
struct FixedInit{T <: Number} <: AbstractInitialPrior
    value::T
end

"""
    GaussianInit(μ0, σ0)

Start each chain with `x(t₁) ~ N(μ0, σ0²)`, independent of the process's stationary
marginal (if any). A large `σ0` leaves the first point nearly free. `σ0` must be
positive, or the constructor throws; for an exact starting value use [`FixedInit`](@ref).
"""
struct GaussianInit{M <: Number, S <: Number} <: AbstractInitialPrior
    μ0::M
    σ0::S
    function GaussianInit{M, S}(μ0, σ0) where {M, S}
        σ0 > 0 || throw(
            ArgumentError("GaussianInit σ0 must be positive; for an exact starting value use FixedInit")
        )
        return new{M, S}(μ0, σ0)
    end
end

GaussianInit(μ0::Number, σ0::Number) = GaussianInit{typeof(μ0), typeof(σ0)}(μ0, σ0)

"""
    UniformInit()

Start each chain uniformly on the circle `(−π, π]`. Valid only for wrapped processes
that are not mean-reverting, such as [`WrappedBrownian`](@ref), whose transitions leave
the uniform law invariant; the chain then absorbs the per-track phase offset. A
mean-reverting circular process ([`WrappedOrnsteinUhlenbeck`](@ref)) uses
[`StationaryInit`](@ref) instead.
"""
struct UniformInit <: AbstractInitialPrior end

"""
    initial_moments(init::AbstractInitialPrior, p::AbstractGaussMarkovProcess) -> (m₀, P₀)

Mean and variance of a chain's first time stamp under `init`.
"""
initial_moments(::StationaryInit, p::AbstractGaussMarkovProcess) = stationary_moments(p)
initial_moments(i::GaussianInit, ::AbstractGaussMarkovProcess) = (i.μ0, abs2(i.σ0))
# Delta limit. A FixedInit first point is always reference-fixed, so the free-point walks
# never use these moments; if they did, the zero variance would give Inf/NaN.
initial_moments(i::FixedInit, ::AbstractGaussMarkovProcess) = (i.value, zero(abs2(i.value)))
initial_moments(::UniformInit, ::AbstractGaussMarkovProcess) = throw(
    ArgumentError("UniformInit has no Gaussian moments; it is only valid for wrapped processes")
)

"""
    marginal_moments(init, p, Δt) -> (m, P)

Moments of the chain marginal a gap `Δt` after the chain start: the initial moments
propagated through the transition law, `m = μ + Φ(m₀ − μ)`, `P = Φ²P₀ + Q`. At `Δt == 0`
this is `initial_moments` exactly; for [`StationaryInit`](@ref) the marginal is the
stationary one at every gap.
"""
function marginal_moments(init::AbstractInitialPrior, p::AbstractGaussMarkovProcess, Δt)
    m₀, P₀ = initial_moments(init, p)
    iszero(Δt) && return (m₀, P₀)
    μ = _process_mean(p)
    Φ, Q = transition_moments(p, Δt)
    return (μ + Φ * (m₀ - μ), Φ^2 * P₀ + Q)
end
marginal_moments(::StationaryInit, p::AbstractGaussMarkovProcess, _) = stationary_moments(p)

# Branchless form for the batched log-density, where `Δt` is traced: `z` is the
# host-computed `iszero(Δt)` indicator.
function marginal_moments(init::AbstractInitialPrior, p::AbstractGaussMarkovProcess, Δt, z)
    m₀, P₀ = initial_moments(init, p)
    μ = _process_mean(p)
    Φ, Q = transition_moments(p, Δt)
    return (z * m₀ + (1 - z) * (μ + Φ * (m₀ - μ)), z * P₀ + (1 - z) * (Φ^2 * P₀ + Q))
end
marginal_moments(::StationaryInit, p::AbstractGaussMarkovProcess, _, _) = stationary_moments(p)

# Center and variance of the Gaussian window for the sheet weight of a wrapped chain's
# first free point. The sheets sum to one for any proper window; the init's own marginal
# makes the flat target reproduce it. `UniformInit` and `FixedInit` (always
# reference-fixed, so unused) have no proper Gaussian marginal and take the unit window at
# zero, which keeps the unused term finite.
@inline _init_window(init::AbstractInitialPrior, p, x1) = initial_moments(init, p)
@inline _init_window(::UniformInit, p, x1) = (zero(x1), one(x1))
@inline _init_window(::FixedInit, p, x1) = (zero(x1), one(x1))

# Log sheet weight of a wrapped chain's first free point over `_init_window`; the weights
# sum to one over the sheets.
@inline function _init_sheet_logweight(init, p, x1)
    m, P = _init_window(init, p, x1)
    δ = x1 - m
    return _normal_logpdf(δ, P) - _wn_logpdf(δ, P)
end

# Validity of an (init, process) pair, checked at site-prior construction. Written against
# the `is_wrapped`/`isstationary` traits, so a new process needs no method here.
# StationaryInit needs a stationary process, GaussianInit a real-line one, and UniformInit
# a nonstationary wrapped one; FixedInit is valid for all.
function _check_init(::StationaryInit, p::AbstractGaussMarkovProcess)
    isstationary(p) || throw(
        ArgumentError(
            "StationaryInit requires a stationary process, but $(nameof(typeof(p))) has " *
                "no stationary distribution. Pick an explicit initial prior, e.g. " *
                (
                is_wrapped(p) ? "UniformInit() or FixedInit(value)." :
                    "GaussianInit(μ0, σ0) or FixedInit(value)."
            )
        )
    )
    return nothing
end
function _check_init(::GaussianInit, p::AbstractGaussMarkovProcess)
    is_wrapped(p) && throw(
        ArgumentError(
            "GaussianInit is not supported for the wrapped process $(nameof(typeof(p))); " *
                "use UniformInit() or FixedInit(value)"
        )
    )
    return nothing
end
function _check_init(::UniformInit, p::AbstractGaussMarkovProcess)
    is_wrapped(p) || throw(
        ArgumentError(
            "UniformInit is only valid for circular (wrapped) processes; " *
                "$(nameof(typeof(p))) lives on the real line"
        )
    )
    isstationary(p) && throw(
        ArgumentError(
            "$(nameof(typeof(p))) reverts to its circular mean, so uniform is not its " *
                "stationary law: use StationaryInit() for the WN(μ, σ²) marginal, or " *
                "FixedInit(value) for an exact start"
        )
    )
    return nothing
end
_check_init(::AbstractInitialPrior, ::AbstractGaussMarkovProcess) = nothing

# Float type the init contributes to the chains' working precision (`_working_type`).
_init_paramtype(::StationaryInit) = Union{}
_init_paramtype(::UniformInit) = Union{}
_init_paramtype(i::FixedInit) = float(typeof(i.value))
_init_paramtype(i::GaussianInit) = float(promote_type(typeof(i.μ0), typeof(i.σ0)))
