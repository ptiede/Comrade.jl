export GaussMarkovSitePrior, siteparams
export Centered, NonCentered, AngleEmbedded

"""
    AbstractChainParameterization

Abstract type for the latent coordinates of a [`GaussMarkovSitePrior`](@ref) chain:
[`Centered`](@ref), [`NonCentered`](@ref), or [`AngleEmbedded`](@ref). See the site
prior's docstring for which to use when.
"""
abstract type AbstractChainParameterization end

"""
    Centered()

Chain coordinates equal to the free points' values: the gain itself or, for a wrapped
process, the raw angle with a sheet weight. The flat target keeps the chain's prior
correlation. Only `asflat` is supported. See [`GaussMarkovSitePrior`](@ref).
"""
struct Centered <: AbstractChainParameterization end

"""
    NonCentered()

Whitened chain coordinates: each free point's coordinate is an iid standard normal
variate, so the flat target has none of the chain's prior correlation. On the real line
the coordinate is mapped to the value through the chain's exact conditional moments,
including the reference-station bridge and the initial marginal. On the circle it is the
step onto the previous phase's transition mean divided by the step's standard deviation,
and the chain's first free point keeps a raw-angle coordinate. Requires a transition mean
that is affine in the previous point or `2π`-equivariant (like the sine drift of
[`VonMisesProcess`](@ref)), which excludes [`WrappedOrnsteinUhlenbeck`](@ref). This is the
default. See [`GaussMarkovSitePrior`](@ref).
"""
struct NonCentered <: AbstractChainParameterization end

"""
    AngleEmbedded()

Chain coordinates that embed each free phase as two latent reals through the angle
transform of `DiagonalVonMises`, doubling the dimension. These coordinates have no `2π`
sheet structure. Valid only for a wrapped process. See [`GaussMarkovSitePrior`](@ref).
"""
struct AngleEmbedded <: AbstractChainParameterization end

"""
    GaussMarkovSitePrior(seg::TimeSegmentation, process::AbstractGaussMarkovProcess; init = StationaryInit(), param = NonCentered())

A site prior correlated in time, following the continuous-time Gauss-Markov `process`
(e.g. [`OrnsteinUhlenbeck`](@ref)) at the times given by the segmentation `seg`. The
log-density factorizes over each site's time chain using the process's exact
discretization, so it costs `O(n)` in the number of times and handles irregular time
sampling exactly. With multiple frequency channels, each (site, channel) pair is an
independent chain.

Process fields that are `Distributions.Distribution`s are fitted hyperparameters. When
any are present the prior sample is a `NamedTuple` `(params = SiteArray, hyperparams =
NamedTuple)`, as for `HierarchicalPrior`; otherwise the sample is a plain `SiteArray`.
Each site gets its own hyperparameters: the process passed to `ArrayPrior` is a per-site
template, and each site's fitted fields are nested under the site's name (e.g.
`x.instrument.lg.hyperparams.AA.σ`). The usual `ArrayPrior` per-site overrides replace the
template, or the whole prior, for that site. `caltable`/`plotcaltable` accept the
hierarchical sample directly; for elementwise postprocessing use the `params` field, e.g.
`exp.(x.instrument.lg.params)`.

# Parameterization

The `param` keyword picks the chain's latent coordinates: [`NonCentered`](@ref),
[`Centered`](@ref), or, for wrapped processes only, [`AngleEmbedded`](@ref).
`centered = true` / `centered = false` is shorthand for `Centered()` / `NonCentered()`;
passing both `param` and `centered` throws. The default is `NonCentered` for every process
that supports it, which is all of them except [`WrappedOrnsteinUhlenbeck`](@ref).

Centered coordinates give the sampler the chain's prior correlation. For a random walk
the condition number of that correlation grows with the number of time stamps, and a
diagonal mass matrix cannot absorb it. A common shift of one station's whole phase chain,
which trades against source position and structure, is then a diagonal direction through
every coordinate of the chain; in non-centered coordinates it is a single coordinate (the
chain's level), and the increments are prior-standardized. `Centered` can mix better when
every point is strongly constrained by the data and the prior is weak.

For a real-line process (`OrnsteinUhlenbeck`, `BrownianMotion`), `NonCentered` also
removes the funnel between fitted hyperparameters and gains and makes the transport to
the Std spaces exact, so both `asflat` and `ascube` work. `Centered` supports `asflat`
only.

For a wrapped (circular) process every parameterization supports `asflat` only. The
chain density is `2π`-periodic, so a real coordinate read directly as an angle gives an
improper flat target, with identical modes spaced `2π` apart in every coordinate.
`Centered` and `NonCentered` therefore multiply each free phase's flat density by the
sheet weight `N(δ; 0, Q) / WN(δ; Q)`, where `δ` is its step from the previous chain value
and `Q` that step's variance. The weights sum to one over the `2π` sheets, so the circular
model is unchanged while the flat target becomes proper; a `2π` step costs
`exp(−(2π)²/2Q)`.

Under `Centered` the flat target is the unwrapped Gaussian random walk in the raw angles.
Under `NonCentered` each free phase's coordinate is its increment divided by the step's
standard deviation, the flat target is `N(0, I)` in those coordinates, and the phase is
the wrapped running sum. In both, the chain's first free point keeps a raw-angle
coordinate. For a non-centered wrapped chain, `logprior` of a constrained sample contains
only the chain-opening terms, and the rest of the density is in the flat transform's
log-Jacobian; only their sum, which samplers and optimizers evaluate, is the chain density.
On the circle `NonCentered` requires a transition mean that carries the previous point's
sheet: [`WrappedBrownian`](@ref) (Φ ≡ 1) or the `2π`-equivariant sine drift of
[`VonMisesProcess`](@ref), whose recurrence is resolved sequentially. The shortest-arc
drift of [`WrappedOrnsteinUhlenbeck`](@ref) does not, so that process is `Centered` only.

[`AngleEmbedded`](@ref) embeds each free phase as two latent reals through the angle
transform of `DiagonalVonMises`. It has no sheet structure (a phase wrap is a continuous
winding of the latent point), at the cost of twice the phase dimension, a per-pair
geometry that a diagonal mass matrix cannot capture, and the chain correlation of the
centered form.

# Phases

For an absolute phase use [`WrappedBrownian`](@ref) with `init = UniformInit()`: the prior
is `2π`-periodic, needs no separate circular offset term, and its coherence-time
posterior is not biased by phase wraps. For a phase held near a level, such as an R/L
gain ratio phase, use [`WrappedOrnsteinUhlenbeck`](@ref), whose `WN(μ, σ²)` marginal keeps
the chain from wandering. A real-line process (e.g. `OrnsteinUhlenbeck` with
`init = FixedInit(0.0)` on top of a circular offset) models the phase as an unwrapped real
and is appropriate only when excursions stay well below `2π`; larger excursions produce
`2π`-shifted posterior modes from the likelihood's periodicity in `exp(iθ)`. The
`phase = true` option of `ArrayPrior` is not supported with this prior and throws.

# Initial distribution

`init` (see [`AbstractInitialPrior`](@ref)) sets the distribution of each chain's first
time stamp. [`StationaryInit`](@ref) (the default) uses the process's stationary
marginal. [`GaussianInit`](@ref) uses an explicit `N(μ0, σ0²)`, e.g. a wide one so the
starting value is nearly free. [`FixedInit`](@ref) conditions the chain to an exact value
at its first time, the continuous-time analogue of `phase = true`. For an unwrapped phase
chain on top of a separate circular offset, `FixedInit(0.0)` removes the level freedom
that would otherwise trade against the offset and, through the likelihood's `2π`
periodicity, produce `2π`-shifted posterior modes. If the reference scheme also fixes a
chain's first point, the two values must agree, or construction throws.

# Reference stations

On the real line a reference-fixed value is a value of the process, and the free points
are drawn from `p(free | fixed)` exactly. Scattered fixed indices, such as those from
`SEFDReference`, are supported.

On the circle a reference-fixed value is a gauge choice. A station's phase is defined
only up to the per-time global phase that the reference scheme removes, so being the
reference says nothing about the station's own gain phase, and conditioning on it would
pull the station's whole chain toward the reference value. A wrapped chain is therefore
built over the time stamps where the station is free, with transitions taken over the
combined gap across skipped stamps: a point after a reference-fixed stamp is an offset
from the last free point. The fixed values still enter the likelihood but carry no prior
term. For [`WrappedBrownian`](@ref) this chain is the exact marginal of the full one.

# Example

```julia
lg ~ ArrayPrior(GaussMarkovSitePrior(IntegSeg(), OrnsteinUhlenbeck(σ = VLBIExponential(0.1), τ = VLBIExponential(2.0))))
gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 3.0, τ = 1.0)); refant = SEFDReference(0.0))
## phase-fluctuation term meant to be combined with a circular per-track offset
dgp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = 1.0); init = FixedInit(0.0)); refant = SEFDReference(0.0))
```

!!! note
    The example uses the `VLBI*` distributions (`VLBIExponential`, `VLBIInverseGamma`, ...)
    for fitted hyperparameters. These accept a scalar traced argument and so work on both
    the CPU and GPU/`Reactant` backends; the plain `Distributions.jl` types work only on
    the CPU.
"""
struct GaussMarkovSitePrior{S <: TimeSegmentation, P <: AbstractGaussMarkovProcess, I <: AbstractInitialPrior, C <: AbstractChainParameterization} <: AbstractSitePrior
    seg::S
    process::P
    init::I
    param::C
end

_asparam(p::AbstractChainParameterization) = p
_asparam(centered::Bool) = centered ? Centered() : NonCentered()

_default_param(p::AbstractGaussMarkovProcess) = has_noncentered_lift(p) ? NonCentered() : Centered()

function GaussMarkovSitePrior(
        seg::TimeSegmentation, process::AbstractGaussMarkovProcess;
        init::AbstractInitialPrior = StationaryInit(),
        param = nothing, centered = nothing
    )
    (param === nothing || centered === nothing) || throw(
        ArgumentError("pass either `param` or the `centered` shorthand, not both")
    )
    sel = param === nothing ? (centered === nothing ? _default_param(process) : centered) : param
    par = _asparam(sel)
    _check_init(init, process)
    _check_param(par, process)
    return GaussMarkovSitePrior(seg, process, init, par)
end

# Written against the process traits, so a new process needs no method here.
_check_param(::AbstractChainParameterization, ::AbstractGaussMarkovProcess) = nothing
function _check_param(::NonCentered, p::AbstractGaussMarkovProcess)
    has_noncentered_lift(p) || throw(
        ArgumentError(
            "$(typeof(p).name.name) does not support the NonCentered parameterization: " *
                "its transition mean takes the shortest arc onto the previous point, which " *
                "drops that point's sheet, so the chain cannot be built as a running " *
                "recurrence on the coordinates. Use param = Centered(), or the smooth-drift " *
                "VonMisesProcess, which supports NonCentered."
        )
    )
    return nothing
end
function _check_param(::AngleEmbedded, p::AbstractGaussMarkovProcess)
    is_wrapped(p) || throw(
        ArgumentError(
            "AngleEmbedded is only valid for a wrapped (circular) process; " *
                "$(typeof(p).name.name) lives on the real line."
        )
    )
    return nothing
end

needs_chain_machinery(::GaussMarkovSitePrior) = true

# ----- per-(site, freq) chain specifications --------------------------------
#
# Built once at `set_array` time and constant during sampling.

struct MarkovChainSpec{P <: AbstractGaussMarkovProcess, I <: AbstractInitialPrior, C <: AbstractChainParameterization, K, V <: AbstractVector{<:Integer}, T <: AbstractVector, F <: NamedTuple, G <: NamedTuple}
    process::P     # template; fitted fields are Distribution placeholders
    init::I        # initial distribution of the chain's first time stamp
    hpsel::Val{K}  # the site whose hyperparameters this chain reads, or nothing if fully fixed
    inds::V        # flat indices into the full parameter vector, ascending in time
    ts::T          # chain times in hours (strictly increasing)
    param::C       # which latent coordinates the chain uses
    free::F        # static per-free-point tables (see `_free_point_tables`)
    fsub::G        # inds/ts restricted to the fixed points, materialized (views do not trace)
end

function MarkovChainSpec(
        process::AbstractGaussMarkovProcess, init::AbstractInitialPrior, hpsel::Val, inds, ts, fixedpos, param::AbstractChainParameterization
    )
    return MarkovChainSpec(
        process, init, hpsel, inds, ts, param,
        _free_point_tables(inds, ts, fixedpos),
        (inds = inds[fixedpos], ts = ts[fixedpos]),
    )
end

# A site with an `IIDSitePrior` override under a GaussMarkov default. Fixed points carry
# no prior term. The flat transport node is cached because rebuilding it per evaluation is
# type-unstable for distributions whose support bounds are runtime values.
struct IIDChainSpec{D, V <: AbstractVector{<:Integer}, N}
    dist::D
    freeinds::V
    fnode::N
end

IIDChainSpec(dist, freeinds) = IIDChainSpec(dist, freeinds, PT.transport_node(dist, PT.TVFlat()))

EnzymeRules.inactive_type(::Type{<:MarkovChainSpec}) = true
EnzymeRules.inactive_type(::Type{<:IIDChainSpec}) = true

struct MarkovGroup{P <: AbstractGaussMarkovProcess, I <: AbstractInitialPrior, R <: AbstractChainParameterization, S, N <: NamedTuple, F <: NamedTuple, C <: NamedTuple}
    process::P     # template of the *first* member; fitted fields are Distribution placeholders
    init::I
    sites::Val{S}  # member sites in group order, or `nothing` when no field is fitted
    param::R
    numvals::N     # per-chain host vectors of the numeric (non-fitted) process fields
    free::F        # concatenated `_free_point_tables`, plus `hpi` (point -> member index)
    chain::C       # concatenated per-time-stamp tables for the chain log-density
    fsub::C        # `chain`, restricted to the reference-fixed time stamps
end

EnzymeRules.inactive_type(::Type{<:MarkovGroup}) = true

# Shared interface of a single chain and a batched group: `process`, `init`, `param`,
# `free`, and `_nfree`.
const ChainUnit = Union{MarkovChainSpec, MarkovGroup}

@inline _selecthp(hp::NamedTuple, ::Val{nothing}) = (;)
@inline _selecthp(hp::NamedTuple, ::Val{K}) where {K} = getproperty(hp, K)

# ----- the chain log-density -------------------------------------------------

# Log-density of the chain's first point under the init marginal propagated a gap Δ0 from
# the chain start (nonzero only for the fixed-subset term of `chain_term`). `z` is the
# host-computed `iszero(Δ0)` indicator, which keeps `marginal_moments` branchless. A
# FixedInit first point is always reference-fixed, so its term is identical in ℓ(all) and
# ℓ(fixed subset) and is set to zero in both.
@inline function _first_point_term(init::AbstractInitialPrior, p, x1, Δ0, z)
    m1, P1 = marginal_moments(init, p, Δ0, z)
    return _marginal_logpdf(p, x1 - m1, P1)
end
@inline _first_point_term(::FixedInit, p, x1, Δ0, z) = zero(abs2(x1))
# The uniform law is invariant under the wrapped transitions, so the constant holds at any Δ0.
@inline _first_point_term(::UniformInit, p, x1, Δ0, z) = -log(oftype(x1, 2π))

function _gm_chain_logpdf(p::AbstractGaussMarkovProcess, init::AbstractInitialPrior, x, inds, ts, Δ0 = zero(eltype(ts)))
    x1 = rgetindex(x, rgetindex(inds, firstindex(inds)))
    ℓ = _first_point_term(init, p, x1, Δ0, oftype(Δ0, iszero(Δ0)))
    μ = _process_mean(p)
    @trace track_numbers = false for i in (firstindex(inds) + 1):lastindex(inds)
        xi = rgetindex(x, rgetindex(inds, i))
        xp = rgetindex(x, rgetindex(inds, i - 1))
        ℓ += _transition_logpdf(p, xi, xp, rgetindex(ts, i) - rgetindex(ts, i - 1), μ)
    end
    return ℓ
end

# ----- fused flat evaluation of non-centered wrapped chains --------------------
#
# For a non-centered wrapped chain, each wrapped-normal transition term of the chain
# density cancels against the `−log WN` sheet weight in `_color_chain_wrapped_noncentered!`,
# so neither is computed. The flat transform returns `log N(zₖ; 0, 1)` per stepping point
# plus the chain-opening init weight, and `chain_term` returns only the chain-opening
# first-point terms.
#
# For these units `logpdf` alone is not the circular chain density; only its sum with the
# flat transform's log-Jacobian is. `asflat` sampling, optimization, and `logdensityof`
# evaluate that sum, but a standalone `logprior` of a constrained sample does not.
_fused_flat(u::ChainUnit) = is_wrapped(u.process) && u.param isa NonCentered

# Concrete unit types so these methods are more specific than the generic ones below.
@inline _fused_chain_term(::UniformInit, u::MarkovChainSpec, x, hp) =
    -_nlive(u) * log(oftype(zero(eltype(x)), 2π))
@inline _fused_chain_term(::UniformInit, u::MarkovGroup, x, hp) =
    -_nlive(u) * log(oftype(zero(eltype(x)), 2π))
@inline _fused_chain_term(::FixedInit, ::MarkovChainSpec, x, hp) = zero(eltype(x))
@inline _fused_chain_term(::FixedInit, ::MarkovGroup, x, hp) = zero(eltype(x))
# `opens` marks each live chain's first stamp. It is host data, so the count folds at
# trace time under Reactant.
@inline _nlive(spec::MarkovChainSpec) = one(eltype(spec.ts))
@inline _nlive(g::MarkovGroup) = sum(g.chain.opens)

# Exact for any init. A wrapped chain conditions on nothing except a FixedInit opening,
# whose term is zero, so there is no fixed-subset subtraction.
@inline function _fused_chain_term(init::AbstractInitialPrior, spec::MarkovChainSpec, x, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    x1 = rgetindex(x, rgetindex(spec.inds, firstindex(spec.inds)))
    Δ0 = zero(eltype(spec.ts))
    return _first_point_term(init, p, x1, Δ0, oftype(Δ0, true))
end
@inline function _fused_chain_term(init::AbstractInitialPrior, g::MarkovGroup, x, hp)
    vals = _group_fieldvals(g, hp)
    tab = g.chain
    ℓ = zero(eltype(x))
    @trace track_numbers = false for i in eachindex(tab.inds)
        p = _point_process(g, vals, rgetindex(tab.hpi, i))
        xi = rgetindex(x, rgetindex(tab.inds, i))
        ℓ += rgetindex(tab.opens, i) *
            _first_point_term(init, p, xi, rgetindex(tab.dt0, i), rgetindex(tab.z0, i))
    end
    return ℓ
end

# log p(free | fixed) = ℓ(all) − ℓ(fixed subset): a Markov process restricted to a subset
# of times is Markov with the composed transitions, and its first marginal is the init
# marginal propagated to the subset's first time (`marginal_moments`). The subtracted term
# depends on the hyperparameters; dropping it would bias their posterior.
@inline function chain_term(spec::MarkovChainSpec, x, hp)
    # A fully reference-fixed chain: the two passes below are identical and cancel.
    length(spec.fsub.inds) == length(spec.inds) && return zero(eltype(x))
    _fused_flat(spec) && return _fused_chain_term(spec.init, spec, x, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    ℓ = _gm_chain_logpdf(p, spec.init, x, spec.inds, spec.ts)
    isempty(spec.fsub.inds) && return ℓ
    return ℓ - _gm_chain_logpdf(p, spec.init, x, spec.fsub.inds, spec.fsub.ts, first(spec.fsub.ts) - first(spec.ts))
end

@inline function chain_term(spec::IIDChainSpec, x, hp)
    # `freeinds` is empty when every point of the site is reference-fixed
    ℓ = zero(eltype(x))
    @trace track_numbers = false for k in eachindex(spec.freeinds)
        ℓ += Dists.logpdf(spec.dist, rgetindex(x, rgetindex(spec.freeinds, k)))
    end
    return ℓ
end

# ----- the full-vector distribution ------------------------------------------

# Distribution over the full instrument parameter vector for a `GaussMarkovSitePrior`:
# independent per-(site, frequency) chains, plus IID terms for `IIDSitePrior` overrides,
# exactly conditioned on the reference-fixed indices. Built by `ObservedArrayPrior`.
struct GaussMarkovChainDist{C <: NamedTuple, G <: Tuple, I <: AbstractVector{<:Integer}, F} <: Dists.ContinuousMultivariateDistribution
    chains::C
    # `chains` batched into `MarkovGroup`s (IID specs unchanged), in spec order. Used for
    # dimension and trait queries and by the traced walkers (`_walk_units`).
    groups::G
    fixedinds::I
    fixedvals::F
    len::Int
    # Concatenated `free.tgt` of the chains resolved by the batched scan, in spec order.
    # Their `(a, c)` coefficients share one pair of vectors and one `affine_scan!`. Each
    # chain's first free point has `a = 0`, so the carry resets at chain boundaries.
    scantgt::Vector{Int}
    # Mask over `scantgt`: 1 for an unwrapped phase to wrap onto `(−π, π]` before
    # scattering, 0 for a real-line value. Stored in the sample element type so the masked
    # wrap does not promote.
    scanwrap::F
    anywrap::Bool
end

# A non-centered chain with a nonlinear mean resolves and scatters its own recurrence
# (`_color_chain_wrapped_recur!`), so it must not claim scan slots.
_scan_tgt(spec::ChainUnit) =
    (spec.param isa NonCentered && has_affine_lift(spec.process)) ? spec.free.tgt : Int[]
_scan_tgt(::IIDChainSpec) = Int[]

_scan_wrap(spec::ChainUnit) = fill(Float64(is_wrapped(spec.process)), length(_scan_tgt(spec)))
_scan_wrap(::IIDChainSpec) = Float64[]

function GaussMarkovChainDist(chains::NamedTuple, fixedinds, fixedvals, len::Int)
    groups = _build_groups(values(chains))
    gs = collect(groups)
    scantgt = reduce(vcat, map(_scan_tgt, gs); init = Int[])
    wrap = reduce(vcat, map(_scan_wrap, gs); init = Float64[])
    scanwrap = convert(typeof(fixedvals), wrap)
    return GaussMarkovChainDist(
        chains, groups, fixedinds, fixedvals, len, scantgt, scanwrap, any(isone, wrap)
    )
end

# The batched units (see `_build_groups`).
chaingroups(d::GaussMarkovChainDist) = d.groups
EnzymeRules.inactive(::typeof(chaingroups), args...) = nothing

# Units the per-evaluation walkers iterate. Groups and chains give the same flat layout and
# values but differ in cost. On an accelerator each chain adds kernel launches, so groups
# are faster. On the host, per-chain loops keep the process loop-invariant, so Enzyme
# accumulates hyperparameter gradients into a scalar; a group gathers parameters per point,
# which costs Enzyme a scatter-add and an allocation per iteration (several times slower).
# The Reactant extension overrides this for traced arrays.
_walk_units(d::GaussMarkovChainDist, x) = values(chainspecs(d))
EnzymeRules.inactive(::typeof(_walk_units), args...) = nothing

Base.length(d::GaussMarkovChainDist) = d.len
# The reference-fixed values are allocated in the working precision (`_working_type`).
Base.eltype(d::GaussMarkovChainDist) = eltype(d.fixedvals)
Dists.sampler(d::GaussMarkovChainDist) = d

chainspecs(d::GaussMarkovChainDist) = d.chains
EnzymeRules.inactive(::typeof(chainspecs), args...) = nothing

struct _ChainFix{X, HP}
    x::X
    hp::HP
end
@inline (hp::_ChainFix)(spec::Union{ChainUnit, IIDChainSpec}) = chain_term(spec, hp.x, hp.hp)

function _chain_logpdf(d::GaussMarkovChainDist, x::AbstractVector, hp::NamedTuple)
    fd = _ChainFix(x, hp)
    ls = map(fd, _walk_units(d, x))
    return sum(ls)
end

Dists.logpdf(d::GaussMarkovChainDist, x::AbstractVector{<:Number}) = _chain_logpdf(d, parent(x), (;))
function Dists._rand!(rng::Random.AbstractRNG, d::GaussMarkovChainDist, x::AbstractVector{<:Real})
    return _rand_chains!(rng, x, d, (;))
end

# ----- the triangular (Rosenblatt) chain map -----------------------------------
#
# The bridge sampler, the whitening and its inverse, and the Std-space transports walk the
# free points with the same conditional: at free point `k`, the Gaussian from the left
# (the transition from the previous value, or the init marginal at a chain opening)
# combined with the next reference-fixed value, which by the Markov property is exact.
# With left prior N(μ + b₁, Q₁) and a transition (Φ₂, Q₂) to the centered fixed value y₂:
#     λ = 1/Q₁ + Φ₂²/Q₂,   mean = μ + (b₁/Q₁ + Φ₂ y₂/Q₂)/λ,   var = 1/λ
# with Φ₂ = 0 when no fixed point follows. The (Φ, Q) form stays regular at Φ = 1, so
# random-walk processes use it too. The walk is branchless over the static `spec.free`
# tables, so under Reactant it compiles to a single while loop.

# Weights of the bridge conditional above: `w` on the left mean shift, `wf` on the fixed
# value, and `s` the conditional standard deviation.
@inline function _bridge_weights(Q₁, Φ₂, Q₂)
    λ = inv(Q₁) + Φ₂^2 / Q₂
    w = inv(Q₁ * λ)
    wf = Φ₂ / (Q₂ * λ)
    s = sqrt(inv(λ))
    return w, wf, s
end

@inline function _bridge_moments(μ, b₁, Q₁, Φ₂, y₂, Q₂)
    w, wf, s = _bridge_weights(Q₁, Φ₂, Q₂)
    return μ + w * b₁ + wf * y₂, s
end

# Conditional moments of free point `k` in affine form: mean `a * y_prevfree + c` and scale
# `s`, where `y_prevfree` is the value of the preceding free point and `c` collects the
# init marginal, the bridge term, and any reference-fixed left neighbor. `a` and `c` depend
# only on hyperparameters, tables, and fixed values, so the coloring is a first-order scan
# (`affine_scan!`). The `hasfix` branch is static per chain. The masked dummy neighbor
# reads require `y`'s dummy slots to be zero (see `_fill_fixed!`).
@inline function _free_moments_affine(p, μ, m₀, P₀, y, free, k)
    Φl, Ql = transition_moments(p, rgetindex(free.dtl, k))
    ml = rgetindex(free.mskl, k)
    Φ₁ = ml * Φl
    Q₁ = ml * Ql + (1 - ml) * P₀
    b₀ = (1 - ml) * (m₀ - μ)
    if free.hasfix
        Φf, Qf = transition_moments(p, rgetindex(free.dtf, k))
        mf = rgetindex(free.mskf, k)
        Φ₂ = mf * Φf
        Q₂ = mf * Qf + (1 - mf)
        y₂ = rgetindex(y, rgetindex(free.fidx, k)) - μ
        w, wf, s = _bridge_weights(Q₁, Φ₂, Q₂)
        γ = Φ₁ * w
        c₀ = μ - γ * μ + w * b₀ + wf * y₂
    else
        γ = Φ₁
        c₀ = μ - γ * μ + b₀
        s = sqrt(Q₁)
    end
    lf = rgetindex(free.lfree, k)
    a = γ * lf
    c = c₀ + γ * (1 - lf) * rgetindex(y, rgetindex(free.lidx, k))
    return a, c, s
end

# Conditional moments with the left value read from a fully realized `y`. Exact for both
# free and reference-fixed left neighbors, since only one of `a` and the masked gather in
# `c` is nonzero.
@inline function _free_moments(p, μ, m₀, P₀, y, free, k)
    a, c, s = _free_moments_affine(p, μ, m₀, P₀, y, free, k)
    return a * rgetindex(y, rgetindex(free.lidx, k)) + c, s
end

# Zero-fill `y` and scatter the reference-fixed values. The zero fill is required: the
# branchless walks read masked dummy positions, and 0 * NaN = NaN. `scatter_values!` keeps
# this one vectorized scatter under Reactant.
function _fill_fixed!(y, d::GaussMarkovChainDist)
    fill!(y, zero(eltype(y)))
    scatter_values!(y, d.fixedinds, d.fixedvals)
    return y
end

# `out[offset + k - 1] = f(k)` for `k in 1:n` in one traced loop. Reactant vectorizes a
# traced loop only when it carries a single value, so each staged array gets its own loop.
function _traced_fill!(f, out, offset, n)
    @trace track_numbers = false for k in 1:n
        rsetindex!(out, f(k), offset + k - 1)
    end
    return out
end

# `sum(f(k) for k in 1:n)` through a staged vector (see `_traced_fill!`).
_traced_sum(f, like, n) = sum(_traced_fill!(f, similar(like, n), 1, n))

# First-order affine scan `c[k] ← a[k] * c[k-1] + c[k]` with `c[0] = 0`, overwriting `c`
# (callers pass discardable buffers; the Reactant overload has the same contract). `a` is
# zero at chain-opening and fixed-left points, which resets the carry with exact gradients.
# ComradeReactantExt overloads this with a log-depth doubling scan: a loop-carried
# recurrence would serialize into a `stablehlo.while`, and chains of 10³–10⁴ points rule
# out an O(n²) closed form.
function affine_scan!(a, c)
    acc = zero(eltype(c))
    for k in eachindex(c)
        acc = muladd(a[k], acc, c[k])
        c[k] = acc
    end
    return c
end

# Static per-free-point tables, in coordinate order: target flat index, left neighbor
# index/gap/mask, and next reference-fixed index/gap/mask. A missing neighbor has mask 0,
# a dummy self index, and a unit gap, and the transition coefficient is multiplied by the
# mask. An `Inf` gap must not be used instead: its τ-derivative is `0 * Inf = NaN`.
function _free_point_tables(inds, ts, fixedpos)
    T = eltype(ts)
    nfix = length(fixedpos)
    tgt = Int[]
    lidx = Int[]
    dtl = T[]
    mskl = T[]
    lfree = T[]
    pofs = Int[]
    fidx = Int[]
    dtf = T[]
    mskf = T[]
    jfix = 1
    for i in eachindex(inds)
        while jfix ≤ nfix && fixedpos[jfix] < i
            jfix += 1
        end
        (jfix ≤ nfix && fixedpos[jfix] == i) && continue
        push!(tgt, inds[i])
        if i > firstindex(inds)
            push!(lidx, inds[i - 1])
            push!(dtl, ts[i] - ts[i - 1])
            push!(mskl, one(T))
            # 1 when the left chain point is free, 0 when it is reference-fixed
            push!(lfree, (jfix > 1 && fixedpos[jfix - 1] == i - 1) ? zero(T) : one(T))
        else
            push!(lidx, inds[i])
            push!(dtl, one(T))
            push!(mskl, zero(T))
            push!(lfree, zero(T))
        end
        # 0-based offset, within this chain's latent segment, of the left neighbor's
        # coordinate when that neighbor is free; otherwise this point's own offset, a
        # static in-bounds index for the masked read. Used by the centered wrapped lift.
        push!(pofs, isone(last(lfree)) ? length(tgt) - 2 : length(tgt) - 1)
        if jfix ≤ nfix
            b = fixedpos[jfix]
            push!(fidx, inds[b])
            push!(dtf, ts[b] - ts[i])
            push!(mskf, one(T))
        else
            push!(fidx, inds[i])
            push!(dtf, one(T))
            push!(mskf, zero(T))
        end
    end
    return (; tgt, lidx, dtl, mskl, lfree, pofs, fidx, dtf, mskf, hasfix = nfix > 0)
end

_nfree(spec::MarkovChainSpec) = length(spec.free.tgt)
_nfree(spec::IIDChainSpec) = length(spec.freeinds)

# ----- batched chain groups ----------------------------------------------------
#
# A `MarkovGroup` concatenates a run of adjacent compatible chains: their
# `_free_point_tables`, chain after chain, plus a per-point index into per-chain process
# parameters. The loop bodies are the per-chain ones run over all points, so under Reactant
# each walker is one set of vectorized ops rather than one per chain. Only adjacent specs
# are grouped, so the flat layout matches the per-chain walkers.
#
# Chains are grouped when they share the process type (which fixes the fitted fields), the
# initial prior (by value), and `param`. Numeric field values may differ, since every field
# is read per point. `hasfix` may differ: with `mskf = 0`, `_bridge_weights(Q₁, 0, 1)` is
# the one-sided conditional.

_nfree(g::MarkovGroup) = length(g.free.tgt)

# The process type fixes which fields are fitted, so `hpsel` need not be in the key.
_groupkey(spec::MarkovChainSpec) = (typeof(spec.process), spec.init, spec.param)

# Field names of `p` that are fitted (Distribution) and that are fixed (numeric).
_fitkeys(p::AbstractGaussMarkovProcess) = keys(hyperprior(p))
_numkeys(p::AbstractGaussMarkovProcess) = filter(!Base.Fix2(in, _fitkeys(p)), propertynames(p))

# Concatenate the members' `_free_point_tables`. `pofs` is a within-chain offset, so it is
# shifted by each chain's base; `hpi` maps each point to its member. A chain's first free
# point has `lfree == 0`, so a shifted `pofs` never crosses a chain boundary.
function _group_free(specs)
    ks = (:tgt, :lidx, :dtl, :mskl, :lfree, :fidx, :dtf, :mskf)
    cols = NamedTuple{ks}(map(k -> reduce(vcat, (getproperty(sp.free, k) for sp in specs)), ks))
    ns = [length(sp.free.tgt) for sp in specs]
    bases = cumsum(ns) .- ns
    pofs = reduce(vcat, (sp.free.pofs .+ b for (sp, b) in zip(specs, bases)))
    hpi = reduce(vcat, (fill(j, n) for (j, n) in enumerate(ns)))
    return (; cols..., pofs, hpi, hasfix = any(sp -> sp.free.hasfix, specs))
end

# Per-time-stamp tables for the batched chain log-density: each point's flat index, its
# predecessor's in the same chain, the gap, and `opens`, which marks chain-opening entries
# whose transition term is replaced by the initial-marginal term. An opening entry reads
# itself as predecessor over a 1 hour dummy gap, so the masked transition stays finite.
# `Δ0` propagates the initial marginal to a subset's first time (see `chain_term`), and
# `z0` is its host-computed `iszero` indicator.
function _group_chain(chains, T)
    inds = Int[]; pinds = Int[]; dts = T[]; opens = T[]; dt0 = T[]; hpi = Int[]
    for (j, ch) in chains
        I, ts, Δ0 = ch
        for i in eachindex(I)
            push!(inds, I[i]); push!(hpi, j)
            if i == firstindex(I)
                push!(pinds, I[i]); push!(dts, one(T)); push!(opens, one(T)); push!(dt0, T(Δ0))
            else
                push!(pinds, I[i - 1]); push!(dts, T(ts[i] - ts[i - 1]))
                push!(opens, zero(T)); push!(dt0, zero(T))
            end
        end
    end
    return (; inds, pinds, dts, opens, dt0, z0 = T.(iszero.(dt0)), hpi)
end

function MarkovGroup(specs::AbstractVector{<:MarkovChainSpec})
    s1 = first(specs)
    T = eltype(s1.free.dtl)
    fitk = _fitkeys(s1.process)
    sites = isempty(fitk) ? nothing : Tuple(_hpsite(sp) for sp in specs)
    # A numeric field shared by every member stays a scalar so arithmetic on it
    # constant-folds; only per-site values become a gathered vector.
    numvals = NamedTuple{_numkeys(s1.process)}(
        map(_numkeys(s1.process)) do k
            vs = [getproperty(sp.process, k) for sp in specs]
            return allequal(vs) ? first(vs) : vs
        end
    )
    # A fully reference-fixed chain's two log-density passes cancel, so it is dropped.
    live = [(j, sp) for (j, sp) in enumerate(specs) if length(sp.fsub.inds) != length(sp.inds)]
    chain = _group_chain(((j, (sp.inds, sp.ts, zero(T))) for (j, sp) in live), T)
    # Only chains with reference-fixed points contribute the conditioning subtraction.
    fsub = _group_chain(
        (
            (j, (sp.fsub.inds, sp.fsub.ts, first(sp.fsub.ts) - first(sp.ts)))
                for (j, sp) in live if !isempty(sp.fsub.inds)
        ),
        T
    )
    return MarkovGroup(
        s1.process, s1.init, Val(sites), s1.param, numvals,
        _group_free(specs), chain, fsub
    )
end

# Read through the field: matching on type-parameter position would break silently if
# the parameters were reordered.
_hpsite(spec::MarkovChainSpec) = _hpkey(spec.hpsel)
_hpkey(::Val{K}) where {K} = K

# Partition the spec tuple into maximal runs of groupable chains, preserving order (and
# hence the flat coordinate layout). `IIDChainSpec`s pass through unchanged.
function _build_groups(specs::Tuple)
    units = Any[]
    run = MarkovChainSpec[]
    flush!() = (isempty(run) || (push!(units, MarkovGroup(run)); run = MarkovChainSpec[]))
    for sp in specs
        if sp isa MarkovChainSpec
            (isempty(run) || isequal(_groupkey(first(run)), _groupkey(sp))) || flush!()
            push!(run, sp)
        else
            flush!()
            push!(units, sp)
        end
    end
    flush!()
    return Tuple(units)
end

# ----- per-point process materialization ---------------------------------------
#
# Fitted fields are stacked into a length-`nchain` vector and gathered per point; fixed
# fields stay host constants. `setproperties` rebuilds the process as `materialize` does,
# so downstream formulas match the per-chain path.

# Must be `vcat(t...)` rather than `reduce(vcat, t)`: the fold concatenates a traced array
# with a traced scalar, which falls back to Base's `typed_vcat` and aborts the trace with
# "Scalar indexing is disallowed" (Reactant ≥ 0.2.280).
@inline _stackvals(t::Tuple) = vcat(t...)

@inline _group_fieldvals(g::MarkovGroup, hp) = _group_fieldvals(g, g.sites, hp)
@inline _group_fieldvals(g::MarkovGroup, ::Val{nothing}, hp) = g.numvals
@inline function _group_fieldvals(g::MarkovGroup, ::Val{S}, hp) where {S}
    fitk = _fitkeys(g.process)
    fit = NamedTuple{fitk}(
        map(k -> _stackvals(map(s -> getproperty(getproperty(hp, s), k), S)), fitk)
    )
    return merge(g.numvals, fit)
end

# Shared fields are scalars (see `MarkovGroup`); per-site fields are gathered.
@inline _gatherval(v::AbstractVector, j) = rgetindex(v, j)
@inline _gatherval(x::Number, _) = x

@inline function _point_process(g::MarkovGroup, vals::NamedTuple{F}, j) where {F}
    return setproperties(g.process, NamedTuple{F}(map(v -> _gatherval(v, j), values(vals))))
end

# A point's process parameters: a single chain materializes once, a group gathers per
# point. `_unit_params` is called once per walker and `_unit_process` once per point.

@inline _unit_params(spec::MarkovChainSpec, hp) = materialize(spec.process, _selecthp(hp, spec.hpsel))
@inline _unit_process(::MarkovChainSpec, p, j) = p
@inline _hpidx(::MarkovChainSpec, free, k) = 1

@inline _unit_params(g::MarkovGroup, hp) = _group_fieldvals(g, hp)
@inline _unit_process(g::MarkovGroup, vals, j) = _point_process(g, vals, j)
@inline _hpidx(::MarkovGroup, free, k) = rgetindex(free.hpi, k)

# ----- the batched chain log-density --------------------------------------------

function _group_chain_logpdf(g::MarkovGroup, tab, x, vals)
    ℓ = zero(eltype(x))
    @trace track_numbers = false for i in eachindex(tab.inds)
        p = _point_process(g, vals, rgetindex(tab.hpi, i))
        μ = _process_mean(p)
        xi = rgetindex(x, rgetindex(tab.inds, i))
        xp = rgetindex(x, rgetindex(tab.pinds, i))
        op = rgetindex(tab.opens, i)
        ℓ += op * _first_point_term(g.init, p, xi, rgetindex(tab.dt0, i), rgetindex(tab.z0, i)) +
            (1 - op) * _transition_logpdf(p, xi, xp, rgetindex(tab.dts, i), μ)
    end
    return ℓ
end

@inline function chain_term(g::MarkovGroup, x, hp)
    isempty(g.chain.inds) && return zero(eltype(x))
    _fused_flat(g) && return _fused_chain_term(g.init, g, x, hp)
    vals = _group_fieldvals(g, hp)
    ℓ = _group_chain_logpdf(g, g.chain, x, vals)
    isempty(g.fsub.inds) && return ℓ
    return ℓ - _group_chain_logpdf(g, g.fsub, x, vals)
end

# ----- exact conditional sampling (OU bridge) -----------------------------------
# Used only by `rand`. `_fill_fixed!` zero-initializes, as the masked reads of
# `_free_moments` require.

struct _ChainSpec{R, X, H}
    rng::R
    x::X
    hp::H
end

(hp::_ChainSpec)(spec) = _rand_chain!(hp.rng, hp.x, spec, hp.hp)

function _rand_chains!(rng::Random.AbstractRNG, x::AbstractVector, d::GaussMarkovChainDist, hp::NamedTuple)
    _fill_fixed!(x, d)
    rc = _ChainSpec(rng, x, hp)
    foreach(rc, chainspecs(d))
    return x
end

function _rand_chain!(rng::Random.AbstractRNG, x, spec::IIDChainSpec, hp)
    @trace track_numbers = false for k in eachindex(spec.freeinds)
        rsetindex!(x, rand(rng, spec.dist), rgetindex(spec.freeinds, k))
    end
    return x
end

function _rand_chain!(rng::Random.AbstractRNG, x, spec::MarkovChainSpec, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    # This Gaussian bridge would draw wrong values for a wrapped process, which needs its
    # own `_rand_chain!` method (as for `WrappedBrownian` below).
    is_wrapped(p) && throw(
        ArgumentError(
            "no chain sampler registered for the wrapped process $(typeof(p).name.name); " *
                "define `_rand_chain!(rng, x, spec::MarkovChainSpec{<:$(typeof(p).name.name)}, hp)`"
        )
    )
    μ = _process_mean(p)
    m₀, P₀ = initial_moments(spec.init, p)
    free = spec.free
    @trace track_numbers = false for k in 1:_nfree(spec)
        m, s = _free_moments(p, μ, m₀, P₀, x, free, k)
        rsetindex!(x, m + s * randn(rng, typeof(m)), rgetindex(free.tgt, k))
    end
    return x
end

# Wrapped Brownian chains are sampled in segments around the fixed points. Before the
# first fixed point is a backward walk (the uniform init makes the time-reversed walk a
# walk), after the last a forward walk, and between consecutive fixed points a Brownian
# bridge whose total 2π winding is drawn from the wrapped-normal mixture.
function _sample_winding(rng::Random.AbstractRNG, Δθ, QT)
    # windings beyond ~4 increment sds contribute negligible mass
    W = max(3, ceil(Int, 4 * sqrt(QT) / (2π)))
    # Shift by the w = 0 exponent (the largest, since |Δθ| ≤ π) so the dominant weight is
    # 1. Unshifted weights all underflow once Δθ²/(2QT) ≳ 745, and the CDF walk would
    # then return -W.
    e0 = abs2(Δθ) / (2QT)
    wts = [exp(e0 - abs2(Δθ + 2π * w) / (2QT)) for w in (-W):W]
    u = rand(rng) * sum(wts)
    c = zero(u)
    for (m, wt) in enumerate(wts)
        c += wt
        u <= c && return m - W - 1
    end
    return W
end

function _rand_chain!(rng::Random.AbstractRNG, x, spec::MarkovChainSpec{<:WrappedBrownian}, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    inds = spec.inds
    ts = spec.ts
    T = float(eltype(x))
    # A site that is the reference at every stamp has an empty chain (`ref_is_gauge`);
    # `_fill_fixed!` has placed its values.
    isempty(inds) && return x

    # positions of the fixed points, in chain order (`fsub.inds` preserves it)
    fx = spec.fsub.inds
    fixpos = Vector{Int}(undef, length(fx))
    j = 1
    for i in eachindex(inds)
        j > length(fx) && break
        if inds[i] == fx[j]
            fixpos[j] = i
            j += 1
        end
    end

    step_walk(prev, Δt) = _wrap_angle(prev + sqrt(_increment_var(p, Δt)) * randn(rng, T))

    if isempty(fixpos)
        x[inds[1]] = T(π) * (2rand(rng, T) - 1)
        for i in 2:length(inds)
            x[inds[i]] = step_walk(x[inds[i - 1]], ts[i] - ts[i - 1])
        end
        return x
    end

    # backward walk before the first fixed point
    for i in (fixpos[1] - 1):-1:1
        x[inds[i]] = step_walk(x[inds[i + 1]], ts[i + 1] - ts[i])
    end
    # bridges between consecutive fixed points, with the winding drawn exactly
    for s in 1:(length(fixpos) - 1)
        fa, fb = fixpos[s], fixpos[s + 1]
        fb - fa == 1 && continue
        ta, tb = ts[fa], ts[fb]
        Δθ = _wrap_angle(x[inds[fb]] - x[inds[fa]])
        w = _sample_winding(rng, Δθ, _increment_var(p, tb - ta))
        e = x[inds[fa]] + Δθ + 2π * w  # unwrapped endpoint of this segment
        cur, tp = x[inds[fa]], ta
        for i in (fa + 1):(fb - 1)
            m = cur + (e - cur) * (ts[i] - tp) / (tb - tp)
            v = _increment_var(p, (ts[i] - tp) * (tb - ts[i]) / (tb - tp))
            cur = m + sqrt(v) * randn(rng, T)
            x[inds[i]] = _wrap_angle(cur)
            tp = ts[i]
        end
    end
    # forward walk after the last fixed point
    for i in (fixpos[end] + 1):length(inds)
        x[inds[i]] = step_walk(x[inds[i - 1]], ts[i] - ts[i - 1])
    end
    return x
end

# Mean-reverting wrapped chains: the Gaussian bridge, with the drift taken through the
# process's `_cond_mean` and every draw wrapped onto `(−π, π]`. The innovation is drawn
# normal rather than wrapped normal and bridge windings are not enumerated; both are valid
# only for `σ ≪ π`. The log-density is exact.
function _rand_chain!(
        rng::Random.AbstractRNG, x,
        spec::MarkovChainSpec{<:Union{WrappedOrnsteinUhlenbeck, VonMisesProcess}}, hp
    )
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ = _process_mean(p)
    m₀, P₀ = initial_moments(spec.init, p)
    free = spec.free
    T = float(eltype(x))
    for k in 1:_nfree(spec)
        Φl, Ql = transition_moments(p, free.dtl[k])
        ml = free.mskl[k]
        # left prior: the transition from the previous value, or the init marginal at a
        # chain opening (`ml == 0`)
        Q₁ = ml * Ql + (1 - ml) * P₀
        b₁ = ml * (_cond_mean(p, x[free.lidx[k]], Φl, μ) - μ) + (1 - ml) * (m₀ - μ)
        if free.hasfix
            Φf, Qf = transition_moments(p, free.dtf[k])
            mf = free.mskf[k]
            Φ₂ = mf * Φf
            Q₂ = mf * Qf + (1 - mf)
            y₂ = mf * _wrap_angle(x[free.fidx[k]] - μ)
            m, s = _bridge_moments(zero(T), b₁, Q₁, Φ₂, y₂, Q₂)
        else
            m, s = b₁, sqrt(Q₁)
        end
        x[free.tgt[k]] = _wrap_angle(μ + m + s * randn(rng, T))
    end
    return x
end

# ----- flat and Std transport nodes: whitening / coloring ----------------------
#
# Non-centered coordinates are iid standard normals mapped (colored) through the triangular
# map above: each free point is `m + s * z`. The map is exact, so it transports the Std
# spaces (`ascube`/`StdNormal`) directly; in flat space the `Σ log s` log-Jacobian is
# added. Centered chains use the values as coordinates and support only the flat space.
# IID override chains use their distribution's own node per free point.

# The angle embedding uses two latent reals per free phase; the others use one.
_flat_dim(spec::ChainUnit) =
    spec.param isa AngleEmbedded ? 2 * _nfree(spec) : _nfree(spec)
_flat_dim(spec::IIDChainSpec) = _nfree(spec) * TV.dimension(PT.transport_node(spec.dist, PT.TVFlat()))
_flat_dim(d::GaussMarkovChainDist) = sum(_flat_dim, chaingroups(d); init = 0)

_std_dim(spec::ChainUnit, space) = _nfree(spec)
_std_dim(spec::IIDChainSpec, space) = _nfree(spec) * PT.dimension(PT.transport_node(spec.dist, space))
_std_dim(d::GaussMarkovChainDist, space) = sum(Base.Fix2(_std_dim, space), chaingroups(d); init = 0)

_std_transportable(spec::ChainUnit) = spec.param isa NonCentered && !is_wrapped(spec.process)
_std_transportable(::IIDChainSpec) = true
function _check_std_transportable(d::GaussMarkovChainDist)
    all(_std_transportable, chaingroups(d)) || throw(
        ArgumentError(
            "This GaussMarkovSitePrior cannot be transported to the Std spaces " *
                "(ascube/StdNormal): Centered and AngleEmbedded coordinates need the " *
                "target log-density, which those transports never evaluate, and a wrapped " *
                "(circular) chain carries a sheet weight that is not a quantile " *
                "transport. Use asflat() (or, for a real-line process, NonCentered())."
        )
    )
    return nothing
end

# --- flat (TransformVariables) chain walkers: thread (logjac, index) over the spec tuple
# by `Base.tail` recursion.

# `index` walks the flat latent vector and `sidx` the shared `(av, cv)` scan buffers (see
# `GaussMarkovChainDist.scantgt`). Staging reads `y` only at reference-fixed positions,
# never another chain's values, so all chains stage first and the caller scans once.
@inline _color_specs_flat!(flag, y, x, index, sidx, av, cv, ::Tuple{}, hp) =
    TV.logjac_zero(flag, eltype(x)), index, sidx
@inline function _color_specs_flat!(flag, y, x, index, sidx, av, cv, specs::Tuple, hp)
    ℓ1, index, sidx = _color_chain_flat!(flag, y, x, index, sidx, av, cv, first(specs), hp)
    ℓ2, index, sidx = _color_specs_flat!(flag, y, x, index, sidx, av, cv, Base.tail(specs), hp)
    return ℓ1 + ℓ2, index, sidx
end

# Shared by the flat and Std colorings: the `(a, c)` buffers for the batched scan, and the
# scan and scatter that resolve them.
_scan_buffers(y, d::GaussMarkovChainDist) =
    (similar(y, length(d.scantgt)), similar(y, length(d.scantgt)))
# Wraps the scanned unwrapped phases of non-centered wrapped chains, leaving real-line
# values alone (see `GaussMarkovChainDist.scanwrap`).
function _wrap_masked(v, m)
    T2π = oftype(float(zero(eltype(v))), 2π)
    return v .- m .* T2π .* round.(v ./ T2π)
end

function _resolve_scan!(y, d::GaussMarkovChainDist, av, cv)
    isempty(d.scantgt) && return y
    v = affine_scan!(av, cv)
    scatter_values!(y, d.scantgt, d.anywrap ? _wrap_masked(v, d.scanwrap) : v)
    return y
end

function _color_all_flat!(flag, y, x, index, d::GaussMarkovChainDist, hp)
    av, cv = _scan_buffers(y, d)
    ℓ, index, _ = _color_specs_flat!(flag, y, x, index, 1, av, cv, _walk_units(d, x), hp)
    _resolve_scan!(y, d, av, cv)
    return ℓ, index
end

# The affine coloring coefficients `(a, c, s)` of free point `k` of a chain unit.
function _unit_affine(spec, vals, y, free, k)
    p = _unit_process(spec, vals, _hpidx(spec, free, k))
    μ = _process_mean(p)
    m₀, P₀ = initial_moments(spec.init, p)
    return _free_moments_affine(p, μ, m₀, P₀, y, free, k)
end

# Angle-embedded wrapped chains: each free phase is two latent reals through the
# `AngleTransform` used by `DiagonalVonMises`. The chain density is evaluated on the angles
# by `chain_term`, and the transform does not depend on hyperparameters.
function _color_chain_wrapped!(flag, y, x, index, spec::ChainUnit)
    free = spec.free
    n = _nfree(spec)
    t = PT.angle_transform()
    # Stage contiguously and scatter once: a loop writing to the scattered `free.tgt`
    # positions serializes under Reactant.
    step = k -> TV.transform_with(flag, t, x, index + 2 * (k - 1))
    θs = _traced_fill!(k -> first(step(k)), similar(y, n), 1, n)
    scatter_values!(y, free.tgt, θs)
    flag isa TV.NoLogJac && return TV.logjac_zero(flag, eltype(x)), index + 2n
    return _traced_sum(k -> step(k)[2], x, n), index + 2n
end

# Centered wrapped chains: one flat coordinate `φ` per free phase, the angle is `φ`
# wrapped to `(−π, π]`, and the flat density includes `_sheet_logweight` so the target is
# proper (see the sheet-weight section of `processes.jl`). The sheet anchor of free point
# `k` is the previous free coordinate (a contiguous read of `x`), the reference-fixed
# value, or at a chain opening the `_init_window`; all are sheet-independent. Reading the
# anchor from `x` rather than `y` keeps this a pointwise loop with no recurrence.
function _color_chain_wrapped_centered!(flag, y, x, index, spec::ChainUnit, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    θs = _traced_fill!(k -> _wrap_angle(rgetindex(x, index + k - 1)), similar(y, n), 1, n)
    scatter_values!(y, free.tgt, θs)
    flag isa TV.NoLogJac && return TV.logjac_zero(flag, eltype(x)), index + n
    ℓ = _traced_sum(x, n) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        φ = rgetindex(x, index + k - 1)
        ml = rgetindex(free.mskl, k)
        lf = rgetindex(free.lfree, k)
        # `free.pofs` is a static index (see `_free_point_tables`).
        prev = rgetindex(x, index + rgetindex(free.pofs, k))
        anchor = lf * prev + (1 - lf) * ml * rgetindex(y, rgetindex(free.lidx, k))
        # Both terms are blended by the mask to keep the loop single-path; the 1 hour
        # dummy gap keeps the unused term finite.
        wt = _sheet_logweight(p, φ, anchor, rgetindex(free.dtl, k), μ)
        w0 = _init_sheet_logweight(spec.init, p, φ)
        return ml * wt + (1 - ml) * w0
    end
    return ℓ, index + n
end

# Non-centered wrapped chains with an affine transition mean (`has_affine_lift`): each
# coordinate is the increment onto the previous chain value divided by the step's standard
# deviation, so the unwrapped phase `φₖ = φₖ₋₁ + √Qₖ zₖ` is resolved by the batched scan.
# The chain-opening point keeps a raw-angle coordinate weighted over `_init_window`. Fused
# with `chain_term` (see `_fused_flat`): each stepping point contributes `log N(zₖ; 0, 1)`.
function _color_chain_wrapped_noncentered!(flag, y, x, index, sidx, av, cv, spec::ChainUnit, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    # `a = lfree` resets the carry wherever the previous chain value is not another free
    # point; the fixed value it steps off is then a constant in `c`.
    _traced_fill!(k -> rgetindex(free.lfree, k), av, sidx, n)
    _traced_fill!(cv, sidx, n) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        _, Q = transition_moments(p, rgetindex(free.dtl, k))
        ml = rgetindex(free.mskl, k)
        lf = rgetindex(free.lfree, k)
        z = rgetindex(x, index + k - 1)
        return ml * (sqrt(Q) * z + (1 - lf) * rgetindex(y, rgetindex(free.lidx, k))) + (1 - ml) * z
    end
    flag isa TV.NoLogJac && return TV.logjac_zero(flag, eltype(x)), index + n, sidx + n
    ℓ = _traced_sum(x, n) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        ml = rgetindex(free.mskl, k)
        z = rgetindex(x, index + k - 1)
        return ml * (-(z^2 + log(oftype(float(z), 2π))) / 2) +
            (1 - ml) * _init_sheet_logweight(spec.init, p, z)
    end
    return ℓ, index + n, sidx + n
end

# Non-centered wrapped chains whose transition mean is `2π`-equivariant but not affine
# (`VonMisesProcess`): the same coordinates and fused target as
# `_color_chain_wrapped_noncentered!`, but the recurrence `θₖ = wrap(m(θₖ₋₁) + √Q zₖ)` is
# nonlinear and is resolved sequentially (one `stablehlo.while` per group under Reactant).
# Equivariance lets the carry be the wrapped angle.
function _color_chain_wrapped_recur!(flag, y, x, index, spec::ChainUnit, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    ℓ = zero(eltype(x))
    θs = similar(y, n)
    θp = zero(eltype(x))
    @trace track_numbers = false for k in 1:n
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        Φ, Q = transition_moments(p, rgetindex(free.dtl, k))
        ml = rgetindex(free.mskl, k)
        lf = rgetindex(free.lfree, k)
        z = rgetindex(x, index + k - 1)
        # The carried free angle, or the reference-fixed value; at a chain opening `ml`
        # masks the step and the raw-angle coordinate is used.
        anchor = lf * θp + (1 - lf) * rgetindex(y, rgetindex(free.lidx, k))
        θk = _wrap_angle(ml * (_cond_mean(p, anchor, Φ, μ) + sqrt(Q) * z) + (1 - ml) * z)
        rsetindex!(θs, θk, k)
        θp = θk
        ℓ += ml * (-(z^2 + log(oftype(float(z), 2π))) / 2) +
            (1 - ml) * _init_sheet_logweight(spec.init, p, z)
    end
    scatter_values!(y, free.tgt, θs)
    flag isa TV.NoLogJac && return TV.logjac_zero(flag, eltype(x)), index + n
    return ℓ, index + n
end

function _color_chain_flat!(flag, y, x, index, sidx, av, cv, spec::ChainUnit, hp)
    free = spec.free
    n = _nfree(spec)
    ℓ0 = TV.logjac_zero(flag, eltype(x))
    n == 0 && return ℓ0, index, sidx
    if spec.param isa AngleEmbedded
        ℓw, iw = _color_chain_wrapped!(flag, y, x, index, spec)
        return ℓw, iw, sidx
    end
    if spec.param isa Centered
        if is_wrapped(spec.process)
            ℓc, ic = _color_chain_wrapped_centered!(flag, y, x, index, spec, hp)
            return ℓc, ic, sidx
        end
        # One vectorized scatter; a scattered-write loop would serialize under Reactant.
        scatter_values!(y, free.tgt, view(x, index:(index + n - 1)))
        return ℓ0, index + n, sidx
    end
    if is_wrapped(spec.process)
        has_affine_lift(spec.process) &&
            return _color_chain_wrapped_noncentered!(flag, y, x, index, sidx, av, cv, spec, hp)
        ℓr, ir = _color_chain_wrapped_recur!(flag, y, x, index, spec, hp)
        return ℓr, ir, sidx
    end
    vals = _unit_params(spec, hp)
    # The recurrence is affine in the previous free value: these loops stage `(a, c)` and
    # the caller's batched `affine_scan!` resolves it.
    coef = k -> _unit_affine(spec, vals, y, free, k)
    _traced_fill!(k -> first(coef(k)), av, sidx, n)
    _traced_fill!(cv, sidx, n) do k
        _, c, s = coef(k)
        return c + s * rgetindex(x, index + k - 1)
    end
    flag isa TV.NoLogJac && return ℓ0, index + n, sidx + n
    return _traced_sum(k -> log(last(coef(k))), x, n), index + n, sidx + n
end

function _color_chain_flat!(flag, y, x, index, sidx, av, cv, spec::IIDChainSpec, hp)
    ℓ = TV.logjac_zero(flag, eltype(x))
    dim = TV.dimension(spec.fnode)
    n = _nfree(spec)
    # Stage contiguously and scatter once (see `_color_chain_wrapped!`).
    step = k -> TV.transform_with(flag, spec.fnode, x, index + (k - 1) * dim)
    ys = _traced_fill!(k -> first(step(k)), similar(y, n), 1, n)
    scatter_values!(y, spec.freeinds, ys)
    flag isa TV.NoLogJac && return ℓ, index + n * dim, sidx
    return _traced_sum(k -> step(k)[2], x, n), index + n * dim, sidx
end

@inline _whiten_specs_flat!(x, index, y, ::Tuple{}, hp) = index
@inline function _whiten_specs_flat!(x, index, y, specs::Tuple, hp)
    index = _whiten_chain_flat!(x, index, y, first(specs), hp)
    return _whiten_specs_flat!(x, index, y, Base.tail(specs), hp)
end

# Inverse of the centered wrapped lift for a `2π`-periodic transition mean
# (`_periodic_mean`): the favored sheet follows from the stored angles alone, so the unwrap
# is pointwise. A chain opening uses the center of `_init_window`.
function _whiten_chain_wrapped_pointwise!(x, index, y, spec::ChainUnit, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    @trace track_numbers = false for k in 1:n
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        m₀, _ = _init_window(spec.init, p, zero(eltype(x)))
        θ = rgetindex(y, rgetindex(free.tgt, k))
        ml = rgetindex(free.mskl, k)
        Φ, _ = transition_moments(p, rgetindex(free.dtl, k))
        yl = rgetindex(y, rgetindex(free.lidx, k))
        m = ml * _cond_mean(p, yl, Φ, μ) + (1 - ml) * m₀
        rsetindex!(x, m + _wrap_angle(θ - m), index + k - 1)
    end
    return index + n
end

# Inverse of the centered wrapped lift for an affine transition mean: each coordinate is
# its sheet anchor plus the wrapped step onto it, which places every point on the sheet the
# weight favors. Anchoring on the previous coordinate makes this a first-order recurrence,
# resolved with `affine_scan!`. `transform ∘ inverse` is the identity on wrapped angles;
# `inverse ∘ transform` is not, since the sheet is not recoverable from the angle.
function _whiten_chain_wrapped_centered!(x, index, y, spec::ChainUnit)
    free = spec.free
    n = _nfree(spec)
    a = similar(x, n)
    c = similar(x, n)
    @trace track_numbers = false for k in 1:n
        θ = rgetindex(y, rgetindex(free.tgt, k))
        ml = rgetindex(free.mskl, k)
        lf = rgetindex(free.lfree, k)
        # zero for a chain-opening point, whose window is centered at zero
        yl = ml * rgetindex(y, rgetindex(free.lidx, k))
        rsetindex!(a, lf, k)
        rsetindex!(c, _wrap_angle(θ - yl) + (1 - lf) * yl, k)
    end
    affine_scan!(a, c)
    @trace track_numbers = false for k in 1:n
        rsetindex!(x, rgetindex(c, k), index + k - 1)
    end
    return index + n
end

# Inverse of the centered wrapped lift for an equivariant nonlinear mean
# (`VonMisesProcess`): `φₖ = m(φₖ₋₁) + wrap(θₖ − m(φₖ₋₁))`. The sheet of the mean depends
# on the previous coordinate, so this is a loop-carried recurrence.
function _whiten_chain_wrapped_centered_recur!(x, index, y, spec::ChainUnit, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    φp = zero(eltype(x))
    @trace track_numbers = false for k in 1:n
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        Φ, _ = transition_moments(p, rgetindex(free.dtl, k))
        m₀, _ = _init_window(spec.init, p, zero(eltype(x)))
        θ = rgetindex(y, rgetindex(free.tgt, k))
        ml = rgetindex(free.mskl, k)
        lf = rgetindex(free.lfree, k)
        # previous coordinate, reference-fixed value, or init window center at an opening
        anchor = lf * φp + (1 - lf) * rgetindex(y, rgetindex(free.lidx, k))
        m = ml * _cond_mean(p, anchor, Φ, μ) + (1 - ml) * m₀
        φk = m + _wrap_angle(θ - m)
        rsetindex!(x, φk, index + k - 1)
        φp = φk
    end
    return index + n
end

# Inverse of the non-centered wrapped lift: each coordinate is the wrapped step onto its
# transition mean divided by the step's standard deviation. Pointwise for every supported
# process, since an equivariant mean modulo 2π depends only on the stored previous angle.
# A chain opening's coordinate is its angle on the sheet nearest its init window's center.
# As for the other wrapped lifts, `inverse ∘ transform` is not the identity.
function _whiten_chain_wrapped_noncentered!(x, index, y, spec::ChainUnit, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    @trace track_numbers = false for k in 1:n
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        Φ, Q = transition_moments(p, rgetindex(free.dtl, k))
        m₀, _ = _init_window(spec.init, p, zero(eltype(x)))
        θ = rgetindex(y, rgetindex(free.tgt, k))
        ml = rgetindex(free.mskl, k)
        m = _cond_mean(p, rgetindex(y, rgetindex(free.lidx, k)), Φ, μ)
        rsetindex!(
            x,
            ml * (_wrap_angle(θ - m) / sqrt(Q)) + (1 - ml) * (m₀ + _wrap_angle(θ - m₀)),
            index + k - 1
        )
    end
    return index + n
end

function _whiten_chain_flat!(x, index, y, spec::ChainUnit, hp)
    free = spec.free
    if spec.param isa AngleEmbedded
        # `(sin θ, cos θ)` with radius 1: transform ∘ inverse is the identity on angles,
        # but inverse ∘ transform does not preserve the radius.
        t = PT.angle_transform()
        dim = TV.dimension(t)
        @trace track_numbers = false for k in 1:_nfree(spec)
            TV.inverse_at!(
                x, index + (k - 1) * dim, t, rgetindex(y, rgetindex(free.tgt, k))
            )
        end
        return index + _nfree(spec) * dim
    end
    if spec.param isa Centered
        n = _nfree(spec)
        if is_wrapped(spec.process)
            # Traits of a concrete process type, so these branches are host-static.
            _periodic_mean(spec.process) &&
                return _whiten_chain_wrapped_pointwise!(x, index, y, spec, hp)
            has_affine_lift(spec.process) &&
                return _whiten_chain_wrapped_centered!(x, index, y, spec)
            return _whiten_chain_wrapped_centered_recur!(x, index, y, spec, hp)
        end
        @trace track_numbers = false for k in 1:n
            rsetindex!(x, rgetindex(y, rgetindex(free.tgt, k)), index + k - 1)
        end
        return index + n
    end
    if is_wrapped(spec.process)
        return _whiten_chain_wrapped_noncentered!(x, index, y, spec, hp)
    end
    vals = _unit_params(spec, hp)
    @trace track_numbers = false for k in 1:_nfree(spec)
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        m₀, P₀ = initial_moments(spec.init, p)
        m, s = _free_moments(p, μ, m₀, P₀, y, free, k)
        rsetindex!(x, (rgetindex(y, rgetindex(free.tgt, k)) - m) / s, index + k - 1)
    end
    return index + _nfree(spec)
end

function _whiten_chain_flat!(x, index, y, spec::IIDChainSpec, hp)
    dim = TV.dimension(spec.fnode)
    @trace track_numbers = false for k in eachindex(spec.freeinds)
        TV.inverse_at!(
            x, index + (k - 1) * dim, spec.fnode,
            rgetindex(y, rgetindex(spec.freeinds, k))
        )
    end
    return index + _nfree(spec) * dim
end

# --- Std-space chain walkers: no log-Jacobian; the latent is converted to a standard
# normal through the space's cdf, as in PT's `ScalarTransport`.

# Same staging/batched-scan structure as the flat walker (see `_color_specs_flat!`).
@inline _color_specs_std!(y, x, index, sidx, av, cv, ::Tuple{}, hp, space) = index, sidx
@inline function _color_specs_std!(y, x, index, sidx, av, cv, specs::Tuple, hp, space)
    index, sidx = _color_chain_std!(y, x, index, sidx, av, cv, first(specs), hp, space)
    return _color_specs_std!(y, x, index, sidx, av, cv, Base.tail(specs), hp, space)
end

function _color_all_std!(y, x, index, d::GaussMarkovChainDist, hp, space)
    av, cv = _scan_buffers(y, d)
    index, _ = _color_specs_std!(y, x, index, 1, av, cv, _walk_units(d, x), hp, space)
    _resolve_scan!(y, d, av, cv)
    return index
end

function _color_chain_std!(y, x, index, sidx, av, cv, spec::ChainUnit, hp, space)
    # centered and wrapped chains are rejected at node construction (`_check_std_transportable`)
    free = spec.free
    n = _nfree(spec)
    n == 0 && return index, sidx
    vals = _unit_params(spec, hp)
    # Same affine split as `_color_chain_flat!`.
    coef = k -> _unit_affine(spec, vals, y, free, k)
    _traced_fill!(k -> first(coef(k)), av, sidx, n)
    _traced_fill!(cv, sidx, n) do k
        _, c, s = coef(k)
        return c + s * _to_stdnormal(space, rgetindex(x, index + k - 1))
    end
    return index + n, sidx + n
end

# Convert a latent coordinate of `space` to a standard normal value and back through the
# space's cdf. For `StdNormal` the map is the identity, avoiding cdf/quantile saturation.
_to_stdnormal(space, z) = PT.space_quantile(PT.StdNormal(), PT._clamp_unit(PT.space_cdf(space, z)))
_to_stdnormal(::PT.StdNormal, z) = z
_from_stdnormal(space, w) = PT.space_quantile(space, PT._clamp_unit(PT.space_cdf(PT.StdNormal(), w)))
_from_stdnormal(::PT.StdNormal, w) = w

function _color_chain_std!(y, x, index, sidx, av, cv, spec::IIDChainSpec, hp, space)
    tin = PT.transport_node(spec.dist, space)
    dim = PT.dimension(tin)
    n = _nfree(spec)
    # Stage contiguously and scatter once (see `_color_chain_wrapped!`).
    ys = similar(y, n)
    @trace track_numbers = false for k in 1:n
        yi, _ = PT.pfwd_step(tin, x, index + (k - 1) * dim)
        rsetindex!(ys, yi, k)
    end
    scatter_values!(y, spec.freeinds, ys)
    return index + n * dim, sidx
end

@inline _whiten_specs_std!(x, index, y, ::Tuple{}, hp, space) = index
@inline function _whiten_specs_std!(x, index, y, specs::Tuple, hp, space)
    index = _whiten_chain_std!(x, index, y, first(specs), hp, space)
    return _whiten_specs_std!(x, index, y, Base.tail(specs), hp, space)
end

function _whiten_chain_std!(x, index, y, spec::ChainUnit, hp, space)
    vals = _unit_params(spec, hp)
    free = spec.free
    n = _nfree(spec)
    @trace track_numbers = false for k in 1:n
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        m₀, P₀ = initial_moments(spec.init, p)
        m, s = _free_moments(p, μ, m₀, P₀, y, free, k)
        rsetindex!(x, _from_stdnormal(space, (rgetindex(y, rgetindex(free.tgt, k)) - m) / s), index + k - 1)
    end
    return index + n
end

function _whiten_chain_std!(x, index, y, spec::IIDChainSpec, hp, space)
    tin = PT.transport_node(spec.dist, space)
    dim = PT.dimension(tin)
    @trace track_numbers = false for k in eachindex(spec.freeinds)
        PT.pback_step!(
            x, index + (k - 1) * dim, tin, rgetindex(y, rgetindex(spec.freeinds, k))
        )
    end
    return index + _nfree(spec) * dim
end

# --- the transport nodes ---

# Flat node for fixed hyperparameters. `ObservedArrayPrior` wraps it in
# `InstrumentTransform`, so it returns the plain full-length vector.
struct MarkovColorTransform{D <: GaussMarkovChainDist} <: TV.VectorTransform
    dists::D
    dim::Int
end

TV.dimension(t::MarkovColorTransform) = t.dim

function TV.transform_with(flag::TV.LogJacFlag, t::MarkovColorTransform, x, index)
    d = t.dists
    y = similar(x, length(d))
    _fill_fixed!(y, d)
    ℓ, index = _color_all_flat!(flag, y, x, index, d, (;))
    return y, ℓ, index
end

function TV.inverse_at!(x::AbstractArray, index, t::MarkovColorTransform, y::AbstractVector)
    return _whiten_specs_flat!(x, index, parent(y), _walk_units(t.dists, parent(y)), (;))
end

TV.inverse_eltype(::MarkovColorTransform, x::Type) = eltype(x)

PT.transport_node(d::GaussMarkovChainDist, ::PT.TVFlat) = MarkovColorTransform(d, _flat_dim(d))

# Std node for the fixed-hyperparameter case; wrapped in `StdInstrumentTransform`.
struct StdMarkovColorTransform{D <: GaussMarkovChainDist, S <: PT.AbstractStdDist} <: PT.AbstractTransport
    dists::D
    space::S
    dim::Int
end

PT.dimension(t::StdMarkovColorTransform) = t.dim

function PT.pfwd_step(t::StdMarkovColorTransform, x, index)
    d = t.dists
    y = similar(x, float(eltype(x)), length(d))
    _fill_fixed!(y, d)
    index = _color_all_std!(y, x, index, d, (;), t.space)
    return y, index
end

function PT.pback_step!(x::AbstractVector, index, t::StdMarkovColorTransform, y)
    return _whiten_specs_std!(x, index, parent(y), _walk_units(t.dists, parent(y)), (;), t.space)
end

function PT.transport_node(d::GaussMarkovChainDist, space::PT.AbstractStdDist)
    _check_std_transportable(d)
    return StdMarkovColorTransform(d, space, _std_dim(d, space))
end

# ----- the hierarchical observed prior ----------------------------------------

# One sample of a hierarchical instrument prior: the gains plus the fitted
# hyperparameters. Methods that dispatch on this shape use this alias.
const HierarchicalSiteSample = NamedTuple{(:params, :hyperparams)}

"""
    siteparams(x)

Extract the instrument parameter values from a single instrument-parameter sample. For
hierarchical samples `(params = SiteArray, hyperparams = NamedTuple)` (produced by
[`GaussMarkovSitePrior`](@ref) with fitted hyperparameters) this returns `x.params`;
for plain `SiteArray` samples it is the identity.
"""
@inline siteparams(x) = x
@inline siteparams(x::HierarchicalSiteSample) = x.params

# Observed prior of an `ArrayPrior` whose `GaussMarkovSitePrior` processes have fitted
# hyperparameters. Samples are `(params = SiteArray, hyperparams = NamedTuple)`, and the
# log-density is the chain density given the hyperparameters plus the hyperprior.
struct ObservedHierarchicalArrayPrior{D <: GaussMarkovChainDist, H <: NamedDist, S <: SiteLookup} <: Dists.ContinuousMultivariateDistribution
    dists::D
    hyperprior::H
    sitemap::S
end

Base.length(d::ObservedHierarchicalArrayPrior) = length(d.dists) + length(d.hyperprior)
Base.eltype(d::ObservedHierarchicalArrayPrior) = eltype(d.dists)

function Dists.logpdf(d::ObservedHierarchicalArrayPrior, x::NamedTuple)
    return _chain_logpdf(d.dists, parent(x.params), x.hyperparams) +
        Dists.logpdf(d.hyperprior, x.hyperparams)
end

function Dists.rand(rng::Random.AbstractRNG, d::ObservedHierarchicalArrayPrior)
    hp = rand(rng, d.hyperprior)
    x = Vector{eltype(d.dists)}(undef, length(d.dists))
    _rand_chains!(rng, x, d.dists, hp)
    return (params = SiteArray(x, d.sitemap), hyperparams = hp)
end

function Dists.rand(rng::Random.AbstractRNG, d::ObservedHierarchicalArrayPrior, n::Int)
    return map(_ -> rand(rng, d), 1:n)
end

# Hyperparameter coordinates come first in the latent vector, so the processes are
# materialized before the chains are mapped. The node's shape depends only on the sitemap,
# the fixed indices, and which fields are fitted, never on hyperparameter values.
struct WhitenedHierarchicalTransform{H, D <: GaussMarkovChainDist, L <: SiteLookup} <: TV.VectorTransform
    hnode::H       # TV node of the hyperprior
    dists::D
    site_map::L
    dim::Int
end

TV.dimension(t::WhitenedHierarchicalTransform) = t.dim

function TV.transform_with(flag::TV.LogJacFlag, t::WhitenedHierarchicalTransform, x, index)
    hp, ℓh, index = TV.transform_with(flag, t.hnode, x, index)
    d = t.dists
    y = similar(x, length(d))
    _fill_fixed!(y, d)
    ℓc, index = _color_all_flat!(flag, y, x, index, d, hp)
    return (params = SiteArray(y, t.site_map), hyperparams = hp), ℓh + ℓc, index
end

function TV.inverse_at!(x::AbstractArray, index, t::WhitenedHierarchicalTransform, y::NamedTuple)
    index = TV.inverse_at!(x, index, t.hnode, y.hyperparams)
    return _whiten_specs_flat!(x, index, parent(y.params), _walk_units(t.dists, parent(y.params)), y.hyperparams)
end

function TV.inverse_eltype(t::WhitenedHierarchicalTransform, ::Type{T}) where {T <: NamedTuple}
    return promote_type(
        eltype(fieldtype(T, :params)), TV.inverse_eltype(t.hnode, fieldtype(T, :hyperparams))
    )
end

function PT.transport_node(d::ObservedHierarchicalArrayPrior, space::PT.TVFlat)
    hnode = PT.transport_node(d.hyperprior, space)
    return WhitenedHierarchicalTransform(
        hnode, d.dists, d.sitemap, TV.dimension(hnode) + _flat_dim(d.dists)
    )
end

struct StdWhitenedHierarchicalTransform{H, D <: GaussMarkovChainDist, L <: SiteLookup, S <: PT.AbstractStdDist} <: PT.AbstractTransport
    hnode::H       # PT node of the hyperprior
    dists::D
    site_map::L
    space::S
    dim::Int
end

PT.dimension(t::StdWhitenedHierarchicalTransform) = t.dim

function PT.pfwd_step(t::StdWhitenedHierarchicalTransform, x, index)
    hp, index = PT.pfwd_step(t.hnode, x, index)
    d = t.dists
    y = similar(x, float(eltype(x)), length(d))
    _fill_fixed!(y, d)
    index = _color_all_std!(y, x, index, d, hp, t.space)
    return (params = SiteArray(y, t.site_map), hyperparams = hp), index
end

function PT.pback_step!(x::AbstractVector, index, t::StdWhitenedHierarchicalTransform, y::NamedTuple)
    index = PT.pback_step!(x, index, t.hnode, y.hyperparams)
    return _whiten_specs_std!(x, index, parent(y.params), _walk_units(t.dists, parent(y.params)), y.hyperparams, t.space)
end

function PT.transport_node(d::ObservedHierarchicalArrayPrior, space::PT.AbstractStdDist)
    _check_std_transportable(d.dists)
    hnode = PT.transport_node(d.hyperprior, space)
    return StdWhitenedHierarchicalTransform(
        hnode, d.dists, d.sitemap, space, PT.dimension(hnode) + _std_dim(d.dists, space)
    )
end

# ----- construction from ArrayPrior --------------------------------------------

# Chain times in absolute UTC hours, the clock of the data's `Ti` field, continued past
# 24 h through `mjd` for multi-day tracks. The origin is the array's first day, shared by
# every chain, so `spec.ts` is comparable to `Ti` and across sites. Only gaps enter the
# density.
function _chain_times(times, inds, mjd0)
    ts = map(i -> (mjd(times[i]) - mjd0) * 24.0 + _center(times[i]), inds)
    for i in (firstindex(ts) + 1):lastindex(ts)
        ts[i] > ts[i - 1] || throw(
            ArgumentError(
                "Chain times for a GaussMarkovSitePrior must be strictly increasing. " *
                    "This should not happen; please file an issue."
            )
        )
    end
    return ts
end

# Split a chain's time stamps into those it is built over and, among them, those it is
# conditioned on. A gauge process (`ref_is_gauge`) drops reference-fixed stamps, so free
# points step over the combined gap. A `FixedInit` opening point is kept and conditioned
# on, since it states where the chain starts.
function _gauge_split(process, init, n::Int, fixedpos)
    ref_is_gauge(process) || return collect(1:n), collect(fixedpos)
    opens = init isa FixedInit && !isempty(fixedpos) && first(fixedpos) == 1
    drop = Set(opens ? fixedpos[2:end] : fixedpos)
    return filter(!in(drop), 1:n), opens ? [1] : Int[]
end

function _chainspec(sp::GaussMarkovSitePrior, site, inds, ts, fixedpos)
    hpk = isempty(hyperprior(sp.process)) ? nothing : site
    keep, kept = _gauge_split(sp.process, sp.init, length(inds), fixedpos)
    return MarkovChainSpec(sp.process, sp.init, Val(hpk), inds[keep], ts[keep], kept, sp.param)
end

function _chainspec(sp::IIDSitePrior, site, inds, ts, fixedpos)
    free = isempty(fixedpos) ? copy(inds) : inds[setdiff(eachindex(inds), fixedpos)]
    return IIDChainSpec(sp.dist, free)
end

# Working float precision of the chains, from the process and IID-override parameters;
# fitted fields contribute their `eltype`.
_site_paramtype(sp::GaussMarkovSitePrior) = promote_type(_paramtype(sp.process), _init_paramtype(sp.init))
_site_paramtype(sp::IIDSitePrior) = float(eltype(sp.dist))
function _working_type(site_dists::NamedTuple)
    T = mapreduce(_site_paramtype, promote_type, values(site_dists); init = Union{})
    return T === Union{} ? Float64 : T
end

function build_markov_observed(d::ArrayPrior, site_dists::NamedTuple, smap::SiteLookup, array)
    d.phase && throw(
        ArgumentError(
            "phase=true is not supported with GaussMarkovSitePrior: the cumulative " *
                "reparameterization it enables is superseded by the correlated prior itself. " *
                "Model phases as unwrapped reals with e.g. an OrnsteinUhlenbeck prior and set " *
                "phase=false."
        )
    )
    d.centroid_station === nothing ||
        throw(ArgumentError("centroid_station is not supported with GaussMarkovSitePrior"))
    # This path has no first-stamp pin for IID sites, so reject one rather than drop it.
    for s in keys(site_dists)
        sd = site_dists[s]
        sd isa IIDSitePrior && initprior(sd) !== nothing && throw(
            ArgumentError(
                "IIDSitePrior init pins are not supported when the same ArrayPrior " *
                    "also carries chain (GaussMarkov) site priors (site $s)"
            )
        )
    end

    T = _working_type(site_dists)
    finds, vals = reference_indices(array, smap, d.refant)
    fixedinds = collect(Int, finds)
    fixedvals = vals === nothing ? T[] : collect(T, vals)

    # FixedInit chains are conditioned to the init value at their first time, in addition
    # to any reference-fixed indices. A reference scheme that also fixes that point must
    # agree with the init value.
    for inds in lookup(smap)
        sp = getproperty(site_dists, smap.sites[first(inds)])
        (sp isa GaussMarkovSitePrior && sp.init isa FixedInit) || continue
        i1 = first(inds)
        v = convert(T, sp.init.value)
        j = findfirst(==(i1), fixedinds)
        if j === nothing
            push!(fixedinds, i1)
            push!(fixedvals, v)
        elseif fixedvals[j] != v
            throw(
                ArgumentError(
                    "Site $(smap.sites[i1]) has FixedInit($(sp.init.value)) but the " *
                        "reference scheme already fixes its chain's first time stamp to " *
                        "$(fixedvals[j]). Make the two values agree or drop one of the " *
                        "constraints."
                )
            )
        end
    end

    # A Vector membership test would make this loop quadratic at IntegSeg scale.
    fixedset = Set(fixedinds)
    # shared origin for every chain (see `_chain_times`)
    mjd0 = minimum(mjd, smap.times)
    chains = map(lookup(smap)) do inds
        s = smap.sites[first(inds)]
        sp = getproperty(site_dists, s)
        ts = _chain_times(smap.times, inds, mjd0)
        fixedpos = findall(in(fixedset), inds)
        return _chainspec(sp, s, collect(Int, inds), ts, fixedpos)
    end

    # Each site's fitted fields, nested under its name. Frequency chains of the same site
    # share its hyperparameters.
    hp = (;)
    for (s, sp) in pairs(site_dists)
        sp isa GaussMarkovSitePrior || continue
        h = hyperprior(sp.process)
        isempty(h) && continue
        hp = merge(hp, NamedTuple{(s,)}((h,)))
    end

    dists = GaussMarkovChainDist(chains, fixedinds, fixedvals, length(smap.times))
    isempty(hp) && return ObservedArrayPrior(dists, smap, false)
    return ObservedHierarchicalArrayPrior(dists, NamedDist(hp), smap)
end
