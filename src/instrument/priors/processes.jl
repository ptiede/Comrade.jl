export AbstractGaussMarkovProcess, OrnsteinUhlenbeck, BrownianMotion, WrappedBrownian
export AbstractInitialPrior, StationaryInit, FixedInit, GaussianInit, UniformInit

"""
    AbstractGaussMarkovProcess

Abstract type for continuous-time, stationary Gauss-Markov processes in state-space form.
These are used to build time-correlated instrument priors (see [`GaussMarkovSitePrior`](@ref))
whose log-density factorizes over the time chain, i.e. is `O(n)` in the number of times.

# Interface

A concrete process `p` must implement

  - `isstationary(p) -> Bool`: whether the process has a stationary distribution.
  - `stationary_moments(p) -> (μ, P)`: the stationary mean and variance (stationary
    processes only).
  - `transition_moments(p, Δt) -> (Φ, Q)`: the exact discretization over a gap `Δt`
    (in hours), so that `x(t + Δt) | x(t) ~ N(μ + Φ*(x(t) - μ), Q)`.

The distribution of each chain's first time stamp is not part of the process: it is set
separately by the [`AbstractInitialPrior`](@ref) of the site prior, which defaults to the
stationary marginal for stationary processes.

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
    isstationary(p::AbstractGaussMarkovProcess) -> Bool

Whether `p` has a stationary distribution, i.e. implements
[`stationary_moments`](@ref Comrade.stationary_moments). [`StationaryInit`](@ref) — the
default chain initialization of [`GaussMarkovSitePrior`](@ref) — requires this. Part of
the process interface; there is no fallback.
"""
function isstationary end

# Whether the process lives on the circle rather than the real line. Wrapped processes
# use density-level transitions (`_transition_logpdf`) instead of Gaussian moments and
# support only the centered parameterization.
is_wrapped(::AbstractGaussMarkovProcess) = false

# Mean entering the transition law `x(t+Δ) | x(t) ~ N(μ + Φ(x(t)−μ), Q)`. Stationary
# processes revert to their stationary mean. Nonstationary processes (Φ ≡ 1) must
# override with any finite value; it cancels out of every Φ = 1 formula.
_process_mean(p::AbstractGaussMarkovProcess) = first(stationary_moments(p))

# Transition log-density over a gap Δt. The default is the Gaussian law implied by
# `transition_moments`; wrapped processes override at the density level.
@inline function _transition_logpdf(p::AbstractGaussMarkovProcess, xi, xp, Δt, μ, T2π)
    Φ, Q = transition_moments(p, Δt)
    return -(abs2(xi - μ - Φ * (xp - μ)) / Q + log(T2π * Q)) / 2
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

# Float element type a materialized sample of `p` lands in: numeric fields contribute their
# own type, fitted (Distribution) fields contribute their `eltype`. Used to derive the
# sample element type of the chain distributions rather than pinning `Float64`.
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
    σ isa Number && σ ≤ 0 && throw(ArgumentError("OrnsteinUhlenbeck σ must be positive"))
    τ isa Number && τ ≤ 0 && throw(ArgumentError("OrnsteinUhlenbeck τ must be positive"))
    μ isa Dists.Distribution && throw(ArgumentError("OrnsteinUhlenbeck μ is not fittable; pass a number"))
    return OrnsteinUhlenbeck(σ, τ, μ)
end

isstationary(::OrnsteinUhlenbeck) = true

@inline stationary_moments(p::OrnsteinUhlenbeck) = (p.μ, p.σ^2)

@inline function transition_moments(p::OrnsteinUhlenbeck, Δt)
    # Keep Φ strictly below 1 so that Q and the conditional (bridge) variances downstream
    # stay positive even when Δt/τ underflows, e.g. during a fitted-τ warmup excursion.
    # Otherwise the chain log-density is NaN instead of merely very negative, and the
    # whitening transform divides by zero.
    Φ = min(exp(-Δt / p.τ), prevfloat(1.0))
    return Φ, p.σ^2 * (1 - Φ^2)
end

"""
    BrownianMotion(; D)

A driftless Brownian motion (Wiener process) with diffusion coefficient `D`: increments
over a gap `Δt` (in hours) are independent `N(0, D*Δt)`, so `√D` is the typical drift per
`√hour`. This is the non-reverting limit of [`OrnsteinUhlenbeck`](@ref) (`τ → ∞` at fixed
diffusion rate `D = 2σ²/τ`), appropriate for quantities that wander without a preferred
level.

The process is nonstationary — the marginal variance grows linearly in time — so a chain
built on it must choose an explicit initial prior ([`GaussianInit`](@ref) or
[`FixedInit`](@ref)); `StationaryInit` is an error.

# Arguments

  - `D`: the diffusion coefficient, in units of the parameter squared per hour. Pass a
    `Number` to fix it or a `Distributions.Distribution` to fit it as a hyperparameter
    with that prior.
"""
struct BrownianMotion{TD} <: AbstractGaussMarkovProcess
    D::TD
end

function BrownianMotion(; D)
    D isa Number && D ≤ 0 && throw(ArgumentError("BrownianMotion D must be positive"))
    return BrownianMotion(D)
end

isstationary(::BrownianMotion) = false
stationary_moments(::BrownianMotion) =
    throw(ArgumentError("BrownianMotion has no stationary distribution"))
_process_mean(p::BrownianMotion) = zero(_paramtype(p))

@inline transition_moments(p::BrownianMotion, Δt) = (one(Δt), p.D * Δt)

"""
    WrappedBrownian(; D)

Brownian motion on the circle: angular increments over a gap `Δt` (hours) are wrapped
normal `WN(0, D*Δt)`, the wrap of [`BrownianMotion`](@ref) onto `(−π, π]`. This is the
natural correlated prior for gain *phases*: the prior is exactly `2π`-periodic, so the
`2π`-shifted posterior modes that a real-line (Gaussian) phase prior produces are exactly
equivalent — which mode a sampler lands in is irrelevant, since the likelihood only sees
`exp(iθ)` — and the diffusion posterior is not biased by phase wraps.

In the visibility domain the process decorrelates exponentially,
`E[e^{i(θ(t+Δ) − θ(t))}] = e^{−DΔ/2}`, so `2/D` is the phase coherence time in hours —
the circular analogue of the damped (Ornstein-Uhlenbeck) correlation. The circular
stationary law is uniform, so with [`UniformInit`](@ref) the chain also absorbs the
per-track phase offset that would otherwise require a separate circular offset term.

A wrapped chain has no Gaussian conditionals, so the whitened parameterization does not
exist. The default is therefore `centered = true`: the raw angles are the flat
coordinates, one per free phase. That lift is kept *proper* — no `2π`-shifted copies of
the posterior, which a raw circular coordinate would otherwise produce — by the sheet
weights described below, which leave the circular model exactly unchanged. Passing
`centered = false` instead embeds each free phase as two latent reals through the same
angle transform used by `DiagonalVonMises`, which has no sheet structure at all but costs
twice the phase dimension; it is the fallback if a fit ever does show sheet-related
pathology. See [`GaussMarkovSitePrior`](@ref) for the trade-off. Only `asflat` is
supported in either form. Reference stations (`refant`) are handled by exact
conditioning, as for the Gaussian processes: wrapped normal transitions compose exactly
(their variances add).

# Arguments

  - `D`: the phase diffusion coefficient in rad² per hour (`2/D` is the coherence time).
    Pass a `Number` to fix it or a `Distributions.Distribution` to fit it as a
    hyperparameter with that prior.

# Example

```julia
gp ~ ArrayPrior(
    GaussMarkovSitePrior(IntegSeg(), WrappedBrownian(D = VLBIExponential(1.0)); init = UniformInit());
    refant = SEFDReference(0.0)
)
```
"""
struct WrappedBrownian{TD} <: AbstractGaussMarkovProcess
    D::TD
end

function WrappedBrownian(; D)
    D isa Number && D ≤ 0 && throw(ArgumentError("WrappedBrownian D must be positive"))
    return WrappedBrownian(D)
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
@inline transition_moments(p::WrappedBrownian, Δt) = (one(Δt), p.D * Δt)

# Wrap an angle to (−π, π] up to floating point; `round` keeps the map traceable and its
# a.e. zero derivative makes the wrap transparent to AD.
@inline _wrap_angle(θ) = θ - oftype(float(θ), 2π) * round(θ / oftype(float(θ), 2π))

# Wrapped-normal log-density of an angular increment Δ with variance Q, using whichever
# theta-function representation converges fast: the image sum over wrapped copies for
# small Q, the Fourier (dual) sum for large Q. Both truncations are exact to machine
# precision at the crossover Q = 4 with three terms. Both branches are evaluated and
# selected with `ifelse` so the expression stays a single trace under Reactant; the dual
# sum's truncation can go slightly negative for small Q, so it is floored before the log
# (the floored branch is never the selected one there).
#
# The image sum is evaluated in log-sum-exp form with its k = 0 term factored out: since
# |Δw| ≤ π, the k = 0 image is always the largest, so every ratio exponent is ≤ 0 and the
# ratio sum lies in [1, 7]. A naive `log(sum(exp(...)))` underflows to `log(0) = -Inf`
# (with a 0/0 = NaN gradient) once abs2(Δw)/(2Q) ≳ 745 — reachable for small-D chains on
# short gaps — where the true log-density is a perfectly representable large negative
# number.
@inline function _wn_logpdf(Δ, Q)
    T2π = oftype(float(Δ), 2π)
    Δw = _wrap_angle(Δ)
    sratio = sum(k -> exp(-(abs2(Δw + T2π * k) - abs2(Δw)) / (2Q)), -3:3)
    limg = -abs2(Δw) / (2Q) + log(sratio) - log(T2π * Q) / 2
    sfouri = ntuple(Val(3)) do j
        Base.@_inline_meta
        exp(-j^2 * Q / 2) * cos(j * Δw)
    end
    sfour = 1 + 2 * sum(sfouri)
    lfour = log(max(sfour, eps(one(sfour)))) - log(T2π)
    return ifelse(Q < 4, limg, lfour)
end

@inline _transition_logpdf(p::WrappedBrownian, xi, xp, Δt, μ, T2π) = _wn_logpdf(xi - xp, p.D * Δt)

# ----- sheet weights: the proper flat lift of a wrapped chain -----------------
#
# A wrapped chain's density is exactly `2π`-periodic in every free angle, so lifting it
# to the line by simply *calling* a real coordinate the angle (the centered
# parameterization) gives an improper flat target: an infinite lattice of identical
# modes spaced `2π` apart in every coordinate. HMC then drifts between sheets, warmup
# variance adaptation blows up, and chains disagree on a quantity the likelihood cannot
# even see.
#
# The lift is made proper by *reweighting the sheets*. A free angle gets the real
# coordinate `φ`, the angle is `φ mod 2π`, and the flat density is multiplied by
#
#     w(φ) = N(φ − ref; 0, Q) / WN(φ − ref; Q)
#
# for a reference `ref` that does not itself depend on the sheet of `φ`. Because the
# wrapped normal is by definition the periodic sum of the normal,
# `Σₘ N(δ + 2πm; 0, Q) = WN(δ; Q)`, the weights of the sheets of any one coordinate sum
# to exactly one:
#
#     Σₘ w(φ + 2πm) = 1.
#
# So the pushforward of the reweighted flat density under `φ ↦ φ mod 2π` is *exactly*
# the chain density — the model is unchanged, and no probability mass is added, moved,
# or lost. What changes is only which representative of each sheet class carries the
# mass: `w` is a smooth partition of unity over the sheets that concentrates on the one
# nearest `ref`.
#
# Taking `ref` to be the previous chain point and `Q` its transition variance makes the
# weight cancel that point's wrapped-normal transition analytically, so the flat target
# is the *unwrapped* Gaussian random walk `Πₖ N(φₖ − φₖ₋₁; 0, Qₖ)` times the (periodic)
# likelihood: proper, smooth, and with a single dominant mode whenever a `2π` step is
# improbable under the process (`exp(−(2π)²/2Q)`).

# Log-density of `N(Δ; 0, Q)`, the unwrapped counterpart of `_wn_logpdf`.
@inline _normal_logpdf(Δ, Q) = -(abs2(Δ) / Q + log(oftype(float(Δ), 2π) * Q)) / 2

"""
    _sheet_logweight(p, xi, xp, Δt, μ, T2π)

Log of the sheet weight of a wrapped chain's free point at `xi`, whose flat lift is
anchored on the (sheet-independent) previous chain value `xp` a gap `Δt` earlier. This is
`log N(Δ; 0, Q) − log WN(Δ; Q)` for the process's own transition law, so summing
`exp` of it over the `2π` sheets of `xi` gives exactly one.
"""
@inline function _sheet_logweight(p::AbstractGaussMarkovProcess, xi, xp, Δt, μ, T2π)
    Φ, Q = transition_moments(p, Δt)
    return _normal_logpdf(xi - μ - Φ * (xp - μ), Q) -
        _transition_logpdf(p, xi, xp, Δt, μ, T2π)
end

# Variance of the Gaussian sheet window placed on a chain-opening free phase, whose
# circular marginal is uniform and so anchors on nothing. It is *not* a prior: the
# `N/WN` ratio sums to one over the sheets for any variance, so the circular law stays
# exactly uniform. The variance only sets how sharply the window picks the principal
# sheet; 1 rad² suppresses the neighbouring sheet by `exp(−(2π)²/2) ≈ 3e-9` while
# staying wide enough to cover `(−π, π]` smoothly.
const _SHEET_WINDOW_VAR = 1.0

"""
    _init_sheet_logweight(init, p, x1)

Log of the sheet weight of a wrapped chain's *first* free point. It has no earlier chain
value to anchor on, so the window is centered at zero; as for `_sheet_logweight`, the
weights of its sheets sum to exactly one, leaving the uniform circular initial law
untouched.
"""
@inline function _init_sheet_logweight(init, p, x1)
    Q = oftype(float(x1), _SHEET_WINDOW_VAR)
    return _normal_logpdf(x1, Q) - _wn_logpdf(x1, Q)
end

# ----- initial distributions p(x(t₁)) ------------------------------------------

"""
    AbstractInitialPrior

Abstract type for the initial distribution `p(x(t₁))` of a [`GaussMarkovSitePrior`](@ref)
chain. The process fixes the transition law between time stamps; the initial prior fixes
how each (site, frequency) chain's *first* time stamp is treated:

  - [`StationaryInit`](@ref): the process's stationary marginal (the default; requires a
    stationary process).
  - [`GaussianInit`](@ref): an explicit normal marginal `N(μ0, σ0²)`, e.g. an
    intentionally wide one.
  - [`FixedInit`](@ref): the first point is held at an exact value by conditioning, so
    the chain describes the evolution relative to a known start.

Initial priors must stay in the Gaussian family (a delta counts as its zero-variance
limit): the exact `O(n)` chain density, the reference-antenna conditioning, and the
whitened transports all rely on closed-form Gaussian conditioning, which a non-Gaussian
first marginal would break.

# Interface

  - `initial_moments(init, p) -> (m₀, P₀)`: mean and variance of the first time stamp.
"""
abstract type AbstractInitialPrior end

"""
    StationaryInit()

Start each chain in the process's stationary marginal `N(μ, P)` from
[`stationary_moments`](@ref Comrade.stationary_moments). This is the natural choice for
a stationary process — every time stamp then has the same marginal — and the default of
[`GaussMarkovSitePrior`](@ref). Requires `isstationary(process)`.
"""
struct StationaryInit <: AbstractInitialPrior end

"""
    FixedInit(value)

Hold each chain exactly at `value` at its first time stamp, so the prior describes the
*evolution* of the quantity relative to a known start. The fixed point is implemented by
exact conditioning of the chain and composes with reference-fixed (`refant`) indices; a
referencing scheme that already fixes a chain's first point to a *different* value is an
error at posterior construction.

`FixedInit(0.0)` on a phase-fluctuation chain modeled on top of a separate circular
offset term removes the level redundancy between chain and offset: without it the
chain's overall level trades off against the offset, and the likelihood's `2π`
periodicity turns that redundancy into spurious posterior modes at `2π`-shifted levels.
"""
struct FixedInit{T <: Number} <: AbstractInitialPrior
    value::T
end

"""
    GaussianInit(μ0, σ0)

Start each chain with `x(t₁) ~ N(μ0, σ0²)`, independent of the process's stationary
marginal (if any). A large `σ0` gives a weakly informative ("diffuse") start in which
the first point is essentially fit freely. `σ0` must be strictly positive; for an exact
starting value use [`FixedInit`](@ref).
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

Start each chain uniformly on the circle `(−π, π]`. Only valid for circular (wrapped)
processes such as [`WrappedBrownian`](@ref), where the uniform law is the process's
stationary circular distribution: it is the natural non-informative start for a phase
and absorbs the per-track phase offset into the chain itself.
"""
struct UniformInit <: AbstractInitialPrior end

"""
    initial_moments(init::AbstractInitialPrior, p::AbstractGaussMarkovProcess) -> (m₀, P₀)

Mean and variance of a chain's first time stamp under `init`.
"""
initial_moments(::StationaryInit, p::AbstractGaussMarkovProcess) = stationary_moments(p)
initial_moments(i::GaussianInit, ::AbstractGaussMarkovProcess) = (i.μ0, abs2(i.σ0))
# The delta limit. A FixedInit chain always has its first point reference-fixed, so the
# free-point walks never consume these moments; if that invariant is ever broken the zero
# variance divides to Inf/NaN rather than silently sampling a wrong marginal.
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

# Validity of an (init, process) pair, checked when the site prior is constructed.
# Written against the `is_wrapped`/`isstationary` traits rather than concrete process
# types, so a new process is fully covered by defining its traits — no per-process
# `_check_init` registrations (which, when forgotten, would accept GaussianInit for a
# wrapped process silently and reject UniformInit with the wrong message):
#   - StationaryInit needs Gaussian stationary moments: stationary and not wrapped.
#   - GaussianInit is a real-line marginal: invalid for wrapped processes.
#   - UniformInit is the circular stationary law: only valid for wrapped processes.
#   - FixedInit (the fallback) is the delta limit and valid everywhere.
function _check_init(::StationaryInit, p::AbstractGaussMarkovProcess)
    is_wrapped(p) && throw(
        ArgumentError(
            "$(nameof(typeof(p))) chains start uniform on the circle: use UniformInit(), " *
                "or FixedInit(value) for an exact start"
        )
    )
    isstationary(p) || throw(
        ArgumentError(
            "StationaryInit requires a stationary process, but $(nameof(typeof(p))) has " *
                "no stationary distribution. Pick an explicit initial prior, e.g. " *
                "GaussianInit(μ0, σ0) or FixedInit(value)."
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
    return nothing
end
_check_init(::AbstractInitialPrior, ::AbstractGaussMarkovProcess) = nothing

# Float type the init contributes to the chains' working precision (`_working_type`).
_init_paramtype(::StationaryInit) = Union{}
_init_paramtype(::UniformInit) = Union{}
_init_paramtype(i::FixedInit) = float(typeof(i.value))
_init_paramtype(i::GaussianInit) = float(promote_type(typeof(i.μ0), typeof(i.σ0)))
