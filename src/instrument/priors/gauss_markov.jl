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

Chain coordinates that are the free points' own values — the gain, or, for a wrapped
process, the raw angle plus its sheet weight. The flat target then carries the chain's
prior correlation, and only `asflat` is supported. See [`GaussMarkovSitePrior`](@ref).
"""
struct Centered <: AbstractChainParameterization end

"""
    NonCentered()

Whitened chain coordinates: each free point's coordinate is an iid standard variate. On
the real line it is colored through the chain's exact conditional moments (including the
reference-station bridge and the initial marginal at the chain opening); on the circle it
is the scaled step onto the previous phase's transition mean, so the phase is the wrapped
running recurrence and the chain-opening point keeps its raw angle. Either way the flat
target carries none of the chain's prior correlation. Requires a transition mean that
carries the previous point's sheet — affine in it, or `2π`-equivariant like the sine
drift of [`VonMisesProcess`](@ref) — which excludes [`WrappedOrnsteinUhlenbeck`](@ref).
The default. See [`GaussMarkovSitePrior`](@ref).
"""
struct NonCentered <: AbstractChainParameterization end

"""
    AngleEmbedded()

Chain coordinates that embed each free phase as two latent reals through the angle
transform of `DiagonalVonMises`, at twice the dimension. The one parameterization of a
wrapped chain with no `2π` sheet structure; valid only for a wrapped process. See
[`GaussMarkovSitePrior`](@ref).
"""
struct AngleEmbedded <: AbstractChainParameterization end

"""
    GaussMarkovSitePrior(seg::TimeSegmentation, process::AbstractGaussMarkovProcess; init = StationaryInit(), param = NonCentered())

A site prior that is correlated in time, following the continuous-time Gauss-Markov
`process` (e.g. [`OrnsteinUhlenbeck`](@ref)) sampled at the times implied by the
segmentation `seg`. The log-density factorizes over each site's time chain using the
process's exact discretization, so the cost is `O(n)` in the number of times and the
irregular time sampling of a VLBI observation is handled exactly. Each site has one chain,
shared by all frequency channels; frequency dependence, such as a delay, goes in the Jones
function through the cell frequency `c.Fr`.

Process fields that are `Distributions.Distribution`s are *fitted* hyperparameters. When
any are present the prior sample is a `NamedTuple` `(params = SiteArray, hyperparams =
NamedTuple)`, mirroring `HierarchicalPrior`; when all fields are numbers the sample is a
plain `SiteArray` like other site priors. Every site gets its *own* hyperparameters: the
process passed to `ArrayPrior` acts as a per-site template, and each site's fitted fields
are nested under the site's name (e.g. `x.instrument.lg.hyperparams.AA.σ`). Per-site
overrides via the usual `ArrayPrior` kwargs change the template (or the prior entirely)
for that site. `caltable`/`plotcaltable` accept the hierarchical sample directly; for
elementwise postprocessing (e.g. `exp.(...)` of log-gains) use the `params` field, i.e.
`exp.(x.instrument.lg.params)`.

The `param` keyword picks the chain's latent coordinates. [`NonCentered`](@ref) whitens
the chain — each free point's coordinate is an iid standard variate, colored through the
exact conditional moments on the real line and, on the circle, the scaled step onto the
previous phase; [`Centered`](@ref) makes the coordinate the value (or angle) itself;
[`AngleEmbedded`](@ref) is a wrapped-only third form. The shorthand `centered = true` / `centered = false` selects
`Centered()` / `NonCentered()`. The default is `NonCentered` for every process that
supports it, which is all of them except [`WrappedOrnsteinUhlenbeck`](@ref).

The default is `NonCentered` because a centered chain hands the sampler the chain's own
prior correlation. For a random walk that correlation is a Green's function whose
condition number grows with the number of time stamps, and a diagonal mass matrix cannot
absorb it. It also matters for what the chain is coupled *to*: a common shift of one
station's whole phase chain trades against source position and structure, and in centered
coordinates that degeneracy is a diagonal ray through every coordinate of the chain. In
non-centered coordinates it is a single coordinate — the chain's level — and the
increments are prior-standardized and nearly orthogonal to the sky. `Centered` can still
mix better when every point is strongly data-constrained and the prior is weak.

For a **real-line process** (`OrnsteinUhlenbeck`, `BrownianMotion`) `NonCentered` colors
iid standard variates through the chain's exact conditional moments. This also removes the
hyperparameter-gain funnel when hyperparameters are fitted and makes the transport to the
Std spaces exact, so both `asflat` and `ascube` work; `Centered` supports `asflat` only.

For a **wrapped (circular) process** every form supports `asflat` only, and the two
non-embedded forms share a *sheet weight*. A wrapped chain's density is exactly
`2π`-periodic, so simply *calling* a real coordinate the angle would make the flat target
improper — an infinite lattice of identical modes, `2π` apart in every coordinate, which a
sampler drifts through and which wrecks warmup variance adaptation. Each free phase is
therefore reweighted by `N(δ; 0, Q) / WN(δ; Q)`, where `δ` is its step onto the previous
chain value and `Q` that step's variance. Those weights sum to exactly one over the `2π`
sheets, so the circular model is untouched — no mass is added, moved, or lost — while the
flat target becomes proper, with a `2π` step costing `exp(−(2π)²/2Q)`.

Under `Centered` the flat target is then the *unwrapped* Gaussian random walk in the raw
angles, correlations and all. Under `NonCentered` the coordinate of each free phase is
instead its increment divided by that step's standard deviation, and the sheet weight
leaves a flat target that is exactly `N(0, I)` in those coordinates; the phase itself is
the wrapped running sum. The chain's first free point keeps a raw-angle coordinate in both
forms — it is the level, and there is nothing before it to take an increment from. The two
halves of a non-centered wrapped chain's flat density are fused for efficiency: `logprior`
of a constrained sample carries only the chain-opening terms, with the remainder in the
flat transform — their sum, which every sampler and optimizer evaluates, is exact and
unchanged.
`NonCentered` requires a transition mean that carries the previous point's sheet, which
on the circle means [`WrappedBrownian`](@ref) (Φ ≡ 1, resolved by the batched affine
scan) or the `2π`-equivariant sine drift of [`VonMisesProcess`](@ref) (nonlinear, so its
recurrence resolves sequentially); the shortest-arc drift of
[`WrappedOrnsteinUhlenbeck`](@ref) wraps its anchor, so that process is `Centered` only.

[`AngleEmbedded`](@ref) embeds each free phase as two latent reals through the same angle
transform used by `DiagonalVonMises`. It is the one parameterization with no sheet
structure at all (a phase wrap is a continuous winding of the latent point), so it is the
safe fallback if a fit ever does show sheet-related pathology — at the cost of twice the
phase dimension, a per-pair rotated-ellipse geometry a diagonal mass matrix cannot
capture, and the chain correlation of the centered form.

This prior is intended for both gain amplitudes and phases. For an *absolute* phase the
preferred process is [`WrappedBrownian`](@ref) with `init = UniformInit()`: the prior is
then exactly `2π`-periodic, needs no separate circular offset term, and its coherence-time
posterior is unbiased by phase wraps. For a phase that is instead pinned near a level —
an R/L gain *ratio* phase, say — use [`WrappedOrnsteinUhlenbeck`](@ref), the circular
mean-reverting process: it is periodic in the same way, but its `WN(μ, σ²)` marginal keeps
the chain from wandering off, which a `WrappedBrownian` chain has nothing to prevent. A
real-line process (e.g. `OrnsteinUhlenbeck` with
`init = FixedInit(0.0)` on top of a circular offset) models the phase as an *unwrapped*
real and is appropriate when the excursions stay well below `2π` — beyond that, the
likelihood's periodicity in `exp(iθ)` produces spurious `2π`-shifted modes that only the
wrapped process removes.

The initial distribution of each chain — how its first time stamp is treated — is set by
`init` (see [`AbstractInitialPrior`](@ref)). [`StationaryInit`](@ref) (the default)
starts the chain in the process's stationary marginal. [`GaussianInit`](@ref) uses an
explicit `N(μ0, σ0²)` first marginal, e.g. an intentionally wide one so the starting
value is essentially fit freely. [`FixedInit`](@ref) conditions the chain to an exact
value at its first time: when the chain models an unwrapped phase fluctuation on top of a
separate (circular) offset term, `FixedInit(0.0)` removes the level freedom that would
otherwise trade off against the offset and, through the likelihood's `2π` periodicity,
produce spurious posterior modes at `2π`-shifted levels. A referencing scheme that
already fixes a chain's first point must agree with the `FixedInit` value; a conflict is
an error.

# Reference stations

How a reference-fixed time stamp enters the chain depends on whether the process lives on
the circle. On the **real line** the fixed value is a value of the process, and the chain
is conditioned on it exactly — the free points are drawn from `p(free | fixed)`, which
works with scattered fixed indices such as those produced by `SensitiveReference`.

On the **circle** it is a gauge choice instead. A station's phase is only defined up to
the per-time global phase the reference scheme removes, so a station being picked as the
reference says nothing about its own gain phase — and conditioning on it would pull the
station's whole chain toward the reference value. A wrapped chain is therefore built over
the time stamps where the station is actually free, with the process's transition law
taken over the composed gap across the skipped stamps: a point following a reference-fixed
stamp is an offset from the last free point, not from the reference value. The fixed
values still enter the likelihood; they just carry no prior term. For
[`WrappedBrownian`](@ref) this free-stamp chain is exactly the marginal of the full one.

# Example

```julia
lg ~ ArrayPrior(GaussMarkovSitePrior(IntegSeg(), OrnsteinUhlenbeck(σ = VLBIExponential(0.1), τ = VLBIExponential(2.0))))
gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 3.0, τ = 1.0)); refant = SensitiveReference(0.0))
## phase-fluctuation term meant to be combined with a circular per-track offset
dgp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = 1.0); init = FixedInit(0.0)); refant = SensitiveReference(0.0))
```

!!! note
    When any hyperparameters are fitted the example above uses the `VLBI*` distributions
    (`VLBIExponential`, `VLBIInverseGamma`, ...). These accept a scalar traced argument and
    so work on both the CPU and GPU/`Reactant` backends; the plain `Distributions.jl` types
    only work on the CPU.
"""
struct GaussMarkovSitePrior{S <: TimeSegmentation, P <: AbstractGaussMarkovProcess, I <: AbstractInitialPrior, C <: AbstractChainParameterization} <: AbstractSitePrior
    seg::S
    process::P
    init::I
    param::C
end

_asparam(p::AbstractChainParameterization) = p
_asparam(centered::Bool) = centered ? Centered() : NonCentered()

# The best coordinates the process supports: increments wherever the chain can be built as
# a running recurrence, values otherwise.
_default_param(p::AbstractGaussMarkovProcess) = has_noncentered_lift(p) ? NonCentered() : Centered()

# `centered` is the two-valued shorthand for `param`; passing both is a conflict rather
# than a precedence rule.
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

# Validity of a (parameterization, process) pair, checked when the site prior is built.
# Written against the `is_wrapped`/`has_affine_lift` traits so a new process is covered by
# defining its traits.
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

# ----- Gauss–Markov chains ------------------------------------------------------
#
# Built once at `set_array` time; everything here is constant during sampling (indices,
# times, and the process *template* whose fitted fields are placeholders).
#
# A `MarkovChains` holds one or more site chains laid out chain after chain:
# the static tables of `_free_point_tables` and `_logpdf_tables`, concatenated, plus a
# per-point index `hpi` into the chains' process parameters. Every walker runs one loop over
# all of its points, so Reactant raises each walker into a single set of vectorized ops
# however many chains it covers. Chains may share a `MarkovChains` when they agree on
# everything the loop bodies dispatch on: the process *type* (which also fixes which fields
# are fitted), the initial prior (by value — `FixedInit`/`GaussianInit` carry numbers the
# moments depend on), and `param`. They need *not* agree on the numeric field values of
# their process templates, which are read per point (see `numvals`). `hasfix` need not
# agree either — a chain with no reference-fixed points carries `mskf = 0`, and
# `_bridge_weights(Q₁, 0, 1)` collapses exactly to the one-sided conditional.
struct MarkovChains{P <: AbstractGaussMarkovProcess, I <: AbstractInitialPrior, R <: AbstractChainParameterization, S, N <: NamedTuple, T <: AbstractVector, F <: NamedTuple, C <: NamedTuple}
    process::P     # template of the first chain; fitted fields are Distribution placeholders
    init::I        # initial distribution of each chain's first time stamp
    sites::Val{S}  # the site whose hyperparameters each chain reads, or nothing if none are fitted
    param::R       # which latent coordinates the chains use
    numvals::N     # numeric process fields: a scalar when every chain shares it, else one per chain
    inds::Vector{Int}               # flat indices into the full parameter vector
    ts::T                           # times in hours, strictly increasing within a chain
    spans::Vector{UnitRange{Int}}   # each chain's positions in `inds` and `ts`
    fixed::Vector{Bool}             # whether each point is reference-fixed
    free::F        # per-free-point tables (see `_free_point_tables`) plus `hpi`
    chain::C       # log-density tables (see `_logpdf_tables`)
    fsub::C        # `chain`, restricted to the reference-fixed points
end

EnzymeRules.inactive_type(::Type{<:MarkovChains}) = true

# One chain over the flat indices `inds` at times `ts`, conditioned on its points at `fixedpos`.
function MarkovChains(
        process::AbstractGaussMarkovProcess, init::AbstractInitialPrior, ::Val{K}, inds, ts, fixedpos, param::AbstractChainParameterization
    ) where {K}
    inds = collect(Int, inds)
    ts = collect(ts)
    fixed = fill(false, length(inds))
    fixed[fixedpos] .= true
    free = _free_point_tables(inds, ts, fixedpos)
    numk = _numkeys(process)
    numvals = NamedTuple{numk}(map(k -> getproperty(process, k), numk))
    # A fully reference-fixed chain conditions on itself: its two log-density passes are
    # identical for any hyperparameters and cancel, so it gets no log-density entries.
    live = length(fixedpos) != length(inds)
    T = eltype(ts)
    chain = _logpdf_tables(live ? [(1, inds, ts, zero(T))] : [], T)
    # Only a chain with reference-fixed points has anything to condition on.
    fsub = _logpdf_tables(
        (live && !isempty(fixedpos)) ? [(1, inds[fixedpos], ts[fixedpos], ts[first(fixedpos)] - first(ts))] : [], T
    )
    return MarkovChains(
        process, init, Val(K === nothing ? nothing : (K,)), param, numvals,
        inds, ts, [firstindex(inds):lastindex(inds)], fixed, (; free..., hpi = fill(1, length(free.tgt))), chain, fsub
    )
end

# Per-chain log-density tables, from `(j, inds, ts, Δ0)` for each chain `j`. `open` has one
# entry per chain: the flat index of its first point, the gap `Δ0` that propagates the
# initial marginal to it (see `chain_term`), and that gap's host-computed `iszero`
# indicator, which keeps `marginal_moments` branchless. `step` has one entry per later
# point: its flat index, its predecessor's, and the gap between them.
function _logpdf_tables(chains, T)
    oinds = Int[]; dt0 = T[]; z0 = T[]; ohpi = Int[]
    sinds = Int[]; pinds = Int[]; dts = T[]; shpi = Int[]
    for (j, I, ts, Δ0) in chains
        push!(oinds, first(I)); push!(dt0, Δ0); push!(z0, iszero(Δ0)); push!(ohpi, j)
        for i in (firstindex(I) + 1):lastindex(I)
            push!(sinds, I[i]); push!(pinds, I[i - 1]); push!(dts, ts[i] - ts[i - 1]); push!(shpi, j)
        end
    end
    return (
        open = (inds = oinds, dt0 = dt0, z0 = z0, hpi = ohpi),
        step = (inds = sinds, pinds = pinds, dts = dts, hpi = shpi),
    )
end

# Sites with an `IIDSitePrior` override mixed under a GaussMarkov default. Fixed points
# carry no prior term, matching `build_dist` which drops them from the product. The flat
# transport node is cached here: rebuilding it per evaluation is type-unstable for
# distributions whose support bounds are runtime values.
struct IIDChainSpec{D, V <: AbstractVector{<:Integer}, N}
    dist::D
    freeinds::V
    fnode::N
end

IIDChainSpec(dist, freeinds) = IIDChainSpec(dist, freeinds, PT.transport_node(dist, PT.TVFlat()))

EnzymeRules.inactive_type(::Type{<:IIDChainSpec}) = true

# ----- the chain log-density -------------------------------------------------

# Log-density of the chain's first point under the init marginal propagated a gap Δ0
# from the chain start (Δ0 = 0 for the full chain; the first fixed time for the
# conditioning subtraction in `chain_term`). `z` is the `iszero(Δ0)` indicator, a static
# table read that keeps `marginal_moments` branchless. For FixedInit the first
# point is always reference-fixed, so its (delta) term appears identically in ℓ(all) and
# ℓ(fixed subset) and is conventionally zero in both.
@inline function _first_point_term(init::AbstractInitialPrior, p, x1, Δ0, z)
    m1, P1 = marginal_moments(init, p, Δ0, z)
    return _marginal_logpdf(p, x1 - m1, P1)
end
@inline _first_point_term(::FixedInit, p, x1, Δ0, z) = zero(abs2(x1))
# Uniform on the circle, which is also invariant under the wrapped transitions, so the
# same constant serves the propagated subset marginal at any Δ0.
@inline _first_point_term(::UniformInit, p, x1, Δ0, z) = -log(oftype(x1, 2π))

# `ℓ(all)` of the chains' points in `tab` (see `_logpdf_tables`), each chain's first point
# under the propagated initial marginal and every later one under its transition law. The
# tables are host data, so the emptiness checks are static; a traced loop over an empty
# range still traces its body, whose table reads would fail.
function _chains_logpdf(u::MarkovChains, tab, x, vals)
    ℓ = _opening_logpdf(u, tab.open, x, vals)
    st = tab.step
    isempty(st.inds) && return ℓ
    @trace track_numbers = false for i in eachindex(st.inds)
        p = _unit_process(u, vals, rgetindex(st.hpi, i))
        xi = rgetindex(x, rgetindex(st.inds, i))
        xp = rgetindex(x, rgetindex(st.pinds, i))
        ℓ += _transition_logpdf(p, xi, xp, rgetindex(st.dts, i), _process_mean(p))
    end
    return ℓ
end

function _opening_logpdf(u::MarkovChains, op, x, vals)
    ℓ = zero(eltype(x))
    isempty(op.inds) && return ℓ
    @trace track_numbers = false for i in eachindex(op.inds)
        p = _unit_process(u, vals, rgetindex(op.hpi, i))
        x1 = rgetindex(x, rgetindex(op.inds, i))
        ℓ += _first_point_term(u.init, p, x1, rgetindex(op.dt0, i), rgetindex(op.z0, i))
    end
    return ℓ
end

# ----- fused flat evaluation of non-centered wrapped chains --------------------
#
# For a non-centered wrapped chain every wrapped-normal transition term of the chain
# density cancels *exactly* against the `−log WN` of the matching sheet weight in
# `_color_chain_wrapped_noncentered!`: both read the same deviation modulo 2π, and WN is
# periodic in it. Evaluating both sides only for them to cancel would cost two
# theta-function sums per free phase in every log-density and gradient evaluation, so the
# boundary between the two halves is moved instead: the coloring returns the whole fused
# flat target — `log N(zₖ; 0, 1)` per stepping point plus the chain-opening init weight —
# and `chain_term` returns only the non-cancelling remainder, the chain-opening
# first-point terms. For the inits a wrapped chain admits that remainder is a per-chain
# constant (`−log 2π` under `UniformInit`, zero under `FixedInit`); the generic method
# keeps any other init exact.
#
# Stated loudly: for these units `logpdf` alone is NOT the circular chain density — only
# its sum with the flat transform's log-density term is meaningful, and that sum is
# unchanged by the fusion. `asflat` sampling, optimization, and `logdensityof` always
# evaluate that sum; a standalone `logprior` of a constrained sample inherits the split.
_fused_flat(u::MarkovChains) = is_wrapped(u.process) && u.param isa NonCentered

# The chain-opening first-point terms, evaluated on the angles. A wrapped gauge chain
# conditions on nothing (`fsub` is empty apart from a FixedInit opening, whose term is zero
# on both sides), so there is no subtraction here. The opening count is host data.
@inline _fused_chain_term(::UniformInit, u::MarkovChains, x, hp) =
    -length(u.chain.open.inds) * log(oftype(zero(eltype(x)), 2π))
@inline _fused_chain_term(::FixedInit, ::MarkovChains, x, hp) = zero(eltype(x))
@inline _fused_chain_term(::AbstractInitialPrior, u::MarkovChains, x, hp) =
    _opening_logpdf(u, u.chain.open, x, _unit_params(u, hp))

# Exact conditioning on the reference-fixed values: the restriction of an order-1 Markov
# process to a subset of times is Markov with the composed transition law over the larger
# gaps, and its first marginal is the init marginal propagated to the subset's first time
# (`marginal_moments`). Hence log p(free | fixed) = ℓ(all) − ℓ(fixed subset). The
# subtracted term depends on the hyperparameters even though the fixed values are
# constants — dropping it would bias the hyperparameter posterior.
@inline function chain_term(u::MarkovChains, x, hp)
    isempty(u.chain.open.inds) && return zero(eltype(x))
    _fused_flat(u) && return _fused_chain_term(u.init, u, x, hp)
    vals = _unit_params(u, hp)
    ℓ = _chains_logpdf(u, u.chain, x, vals)
    isempty(u.fsub.open.inds) && return ℓ
    return ℓ - _chains_logpdf(u, u.fsub, x, vals)
end

@inline function chain_term(spec::IIDChainSpec, x, hp)
    # `init`: `freeinds` is empty when every point of the site is reference-fixed
    ℓ = zero(eltype(x))
    @trace track_numbers = false for k in eachindex(spec.freeinds)
        ℓ += Dists.logpdf(spec.dist, rgetindex(x, rgetindex(spec.freeinds, k)))
    end
    return ℓ
end

# ----- the full-vector distribution ------------------------------------------

#    GaussMarkovChainDist
#
#The observed distribution over the full instrument parameter vector implied by
#[`GaussMarkovSitePrior`](@ref): a set of independent per-site Gauss-Markov
#chains (plus IID terms for `IIDSitePrior` overrides), exactly conditioned on the
#reference-fixed indices. This is an internal type constructed by `ObservedArrayPrior`.
struct GaussMarkovChainDist{C <: NamedTuple, G <: Tuple, I <: AbstractVector{<:Integer}, F} <: Dists.ContinuousMultivariateDistribution
    chains::C
    # Runs of compatible `chains` concatenated into one `MarkovChains` each (plus any
    # `IIDChainSpec`s, unchanged), in chain order, used for the dimension/trait queries and
    # by the traced walkers (see `_walk_units`). `chains`, one `MarkovChains` per chain, is
    # what the host walkers and the sampler iterate.
    groups::G
    fixedinds::I
    fixedvals::F
    len::Int
    # Concatenated `free.tgt` of the non-centered (scanning) Markov specs, in spec
    # order: the flat/std coloring stages every such chain's `(a, c)` coefficients
    # into one shared pair of vectors, resolves all recurrences with a SINGLE
    # `affine_scan!`, and scatters once. Exact because each chain's first free point
    # has `a = 0` (chain-opening or fixed-left), so the carry resets at every chain
    # boundary of the concatenation.
    scantgt::Vector{Int}
    # Elementwise mask over `scantgt`: 1 where the scanned value is an unwrapped phase
    # that must be wrapped onto `(−π, π]` before it is scattered, 0 where it is a value on
    # the real line that must be left alone. Held in the sample element type so the masked
    # wrap is a single promotion-free broadcast.
    scanwrap::F
    anywrap::Bool
end

# Only affine-lift chains participate in the batched scan; a non-centered chain with a
# nonlinear (equivariant) mean resolves its own recurrence and scatters its angles itself
# (`_color_chain_wrapped_recur!`), so it must not claim scan slots it never stages.
_scan_tgt(spec::MarkovChains) =
    (spec.param isa NonCentered && has_affine_lift(spec.process)) ? spec.free.tgt : Int[]
_scan_tgt(::IIDChainSpec) = Int[]

_scan_wrap(spec::MarkovChains) = fill(Float64(is_wrapped(spec.process)), length(_scan_tgt(spec)))
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

# The concatenated units (see `_build_groups`).
chaingroups(d::GaussMarkovChainDist) = d.groups
EnzymeRules.inactive(::typeof(chaingroups), args...) = nothing

# Which units the per-evaluation walkers iterate. A group is the concatenation of its
# chains *in order*, so both choices produce the same flat layout, the same `scantgt`,
# and the same numbers — they differ only in cost model, and the two cost models disagree:
#
#   - On an accelerator every chain costs a kernel launch, and the arithmetic is free at
#     these sizes. Batching makes the launch count O(1) in the number of chains instead of
#     O(nchains) (measured post-XLA: 670 -> 124 kernels at 7 sites, and flat rather than
#     linear as sites are added).
#   - On the host the loops are plain Julia and launches do not exist. A single chain
#     keeps the process loop-invariant (see `_unit_params`), so `σ²`/`stationary_moments`
#     hoist out of the loop and Enzyme accumulates `∂/∂σ` into a *scalar* shadow. A group
#     of several sites reads each point's parameters through a gather, which makes that a
#     scatter-add into a length-nchain shadow and costs Enzyme a per-iteration allocation:
#     measured 3.4x slower gradients (34 -> 118 μs, 21 -> 76 allocs) for a
#     fitted-hyperparameter OU chain.
#
# So: groups where launches dominate, chains where hoisting does. The Reactant extension
# overrides this for traced arrays; everything else takes the host path.
_walk_units(d::GaussMarkovChainDist, x) = values(chainspecs(d))
EnzymeRules.inactive(::typeof(_walk_units), args...) = nothing

Base.length(d::GaussMarkovChainDist) = d.len
# Sample element type follows the reference-fixed values, which `build_markov_observed`
# allocates in the working precision derived from the process/override parameters (see
# `_working_type`). This is `Float64` for the usual `Float64` inputs, but is no longer
# pinned there.
Base.eltype(d::GaussMarkovChainDist) = eltype(d.fixedvals)
Dists.sampler(d::GaussMarkovChainDist) = d

chainspecs(d::GaussMarkovChainDist) = d.chains
EnzymeRules.inactive(::typeof(chainspecs), args...) = nothing

struct _ChainFix{X, HP}
    x::X
    hp::HP
end
@inline (hp::_ChainFix)(spec::Union{MarkovChains, IIDChainSpec}) = chain_term(spec, hp.x, hp.hp)

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
# Every per-point operation on a chain — the bridge sampler, the whitened coloring, its
# inverse, and the Std-space transports — walks the *free points only* and needs the same
# exact conditional moments: at free point `k`, condition on a Gaussian prior-from-the-
# left (the transition from the previously realized value, or the init marginal when the
# point opens its chain) and the *next reference-fixed* value, which by the Markov
# property is the exact conditional. In precision form, with prior N(μ + b₁, Q₁) and a
# fixed-point transition (Φ₂, Q₂) to the centered value y₂:
#     λ = 1/Q₁ + Φ₂²/Q₂,   mean = μ + (b₁/Q₁ + Φ₂ y₂/Q₂)/λ,   var = 1/λ
# with Φ₂ = 0 when no fixed point follows, recovering the one-sided conditionals. The
# (Φ, Q) form stays regular at Φ = 1, so nonstationary (random-walk) processes share the
# loop; the (μ, P, Φ)-only stationary form is singular there. The walk is branchless over
# the static `spec.free` tables so the same loop serves the CPU and, via `@trace`,
# compiles to a single while loop under Reactant (the fixed-point skip logic of a
# sequential walk cannot be loop-carried in a traced while).

# Precision-form weights of the two-sided bridge conditional, N(μ + b₁, Q₁) from the
# left and a fixed-point transition (Φ₂, Q₂) to the centered value y₂ on the right:
# λ = 1/Q₁ + Φ₂²/Q₂, with `w` the weight on the left mean shift, `wf` the weight on the
# fixed value, and `s` the conditional scale. The (Φ, Q) form stays regular at Φ = 1
# and at extreme correlation times (see the transition_moments clamping tests).
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

# Conditional moments of free point `k` in AFFINE form: `m = a * y_prevfree + c` with
# scale `s`, where `y_prevfree` is the realized value of the immediately preceding FREE
# point (the loop-carried value of the coloring recurrence) and `c` collects everything
# else — the init marginal, the next-reference-fixed bridge term, and left neighbors
# that are reference-fixed constants (their masked gather from `y` reads values placed
# by `_fill_fixed!`). Splitting on the static `free.lfree` mask is what lets the
# coloring run as a vectorized first-order scan (`affine_scan!`) instead of a serialized
# recurrence: `a` and `c` depend only on hyperparameters, tables, and fixed values.
# The `hasfix` branch is static (per chain), so chains without reference-fixed points
# skip the Φ₂ transition entirely — halving the transcendental cost for refant-free
# priors. Chain-opening points (left mask 0) take the init marginal (m₀, P₀) as their
# prior-from-the-left; the dummy neighbor reads are scaled by exact static zero masks,
# which requires `y`'s dummy slots to be zero-initialized (see `_fill_fixed!`).
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

# Conditional moments with the left value read from a fully realized `y` (whitening,
# sampling, and the std walkers): `a * y[lidx] + c` is exact for both a free and a
# reference-fixed left neighbor, since exactly one of `a` and the masked gather in `c`
# is nonzero.
@inline function _free_moments(p, μ, m₀, P₀, y, free, k)
    a, c, s = _free_moments_affine(p, μ, m₀, P₀, y, free, k)
    return a * rgetindex(y, rgetindex(free.lidx, k)) + c, s
end

# Zero-initialize the full parameter vector and scatter the reference-fixed values. The
# zero fill matters: the branchless walk reads masked dummy neighbor positions, which
# must not be uninitialized memory (0 * NaN = NaN). `scatter_values!` keeps this a single
# vectorized `stablehlo.scatter` under Reactant (the indices are static and unique) —
# this runs on every logdensity/gradient evaluation, and an element-wise traced loop
# would serialize into a stablehlo.while there.
function _fill_fixed!(y, d::GaussMarkovChainDist)
    fill!(y, zero(eltype(y)))
    scatter_values!(y, d.fixedinds, d.fixedvals)
    return y
end

# `outs[j][offsets[j] + k - 1] = f(k)[j]` for `k in 1:n`, in one traced loop. Reactant's
# loop raising vectorizes it because each output is written once per iteration at the loop
# index. A value carried between iterations would keep it a `stablehlo.while`, so per-point
# log terms are staged as an output and summed afterward.
function _traced_fill!(f, outs::Tuple, offsets::Tuple, n)
    @trace track_numbers = false for k in 1:n
        vals = f(k)
        map((o, v, i) -> rsetindex!(o, v, i + k - 1), outs, vals, offsets)
    end
    return outs
end

_traced_fill!(f, out, offset, n) = only(_traced_fill!(k -> (f(k),), (out,), (offset,), n))

_logjac_sum(flag::TV.NoLogJac, ℓs) = TV.logjac_zero(flag, eltype(ℓs))
_logjac_sum(::TV.LogJac, ℓs) = sum(ℓs)

# First-order affine scan `c[k] ← a[k] * c[k-1] + c[k]` with implicit `c[0] = 0`,
# overwriting `c` (both call sites pass discardable staging temps; the Reactant overload
# shares this contract). The coloring coefficients are exactly zero at chain-opening and
# fixed-left points (`a = γ·lfree` with static zero masks), so the carry resets at segment
# boundaries with no special casing — and since the recurrence composes with plain
# mul/add, those zeros also have exact gradients. Sequential on the CPU;
# ComradeReactantExt overloads this with a log-depth doubling scan, because Reactant's
# loop raising never vectorizes a loop-carried recurrence (it would serialize into a
# `stablehlo.while` in every logdensity/gradient evaluation) and IntegSeg-scale chains
# (10³-10⁴ points) rule out any O(n²) closed form.
function affine_scan!(a, c)
    acc = zero(eltype(c))
    for k in eachindex(c)
        acc = muladd(a[k], acc, c[k])
        c[k] = acc
    end
    return c
end

# Static per-free-point tables for a Markov chain, computed once at construction: for
# free point `k` (in coordinate order) the target flat index, the left-neighbor flat
# index/gap/mask, and the next-reference-fixed flat index/gap/mask. A missing neighbor
# is a `0.0` mask with a dummy self index and unit gap: the transition coefficient is
# *multiplied by the mask*, so the dummy read is scaled by an exact static zero. (An
# `Inf` gap would also zero the primal, but its τ-derivative is `0 * Inf = NaN` under
# AD.)
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
            # 1 when the left chain point is itself free (i.e. the previously colored
            # free point), 0 when it is reference-fixed: the affine-scan coloring needs
            # to know statically which left reads are recurrence reads vs constants.
            push!(lfree, (jfix > 1 && fixedpos[jfix - 1] == i - 1) ? zero(T) : one(T))
        else
            push!(lidx, inds[i])
            push!(dtl, one(T))
            push!(mskl, zero(T))
            push!(lfree, zero(T))
        end
        # Offset, within this chain's contiguous latent segment, of the free point whose
        # coordinate is this one's *left* neighbour. Free points take one coordinate each
        # in coordinate order, so that is `k - 2` (0-based) when the left chain point is
        # itself free. When it is not (reference-fixed, or the chain opening), the read is
        # masked off, so the offset points at this point's own coordinate: a static,
        # always-in-bounds index, which a computed `max(k, 2)` would not be. Only the
        # centered wrapped lift uses this.
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

_nfree(u::MarkovChains) = length(u.free.tgt)
_nfree(spec::IIDChainSpec) = length(spec.freeinds)

# ----- concatenating chains ------------------------------------------------------

# Chains share a `MarkovChains` iff they agree on what the loop bodies dispatch on. `init`
# and `param` compare by value; the process compares by *type*, which also fixes which
# fields are `Distribution`s and hence whether `sites` is `nothing`.
_groupkey(u::MarkovChains) = (typeof(u.process), u.init, u.param)

# Field names of `p` that are fitted (Distribution) and that are fixed (numeric).
_fitkeys(p::AbstractGaussMarkovProcess) = keys(hyperprior(p))
_numkeys(p::AbstractGaussMarkovProcess) = filter(!Base.Fix2(in, _fitkeys(p)), propertynames(p))

_nchains(u::MarkovChains) = length(u.spans)

# Concatenate the columns of per-unit tables, adding `shifts[k][j]` to column `k` of unit `j`.
function _cattable(tabs, shifts::NamedTuple)
    ks = keys(first(tabs))
    return NamedTuple{ks}(
        map(ks) do k
            cols = [getproperty(t, k) for t in tabs]
            haskey(shifts, k) || return reduce(vcat, cols)
            return reduce(vcat, map((c, b) -> c .+ b, cols, shifts[k]))
        end
    )
end

_bases(ns) = cumsum(ns) .- ns

# `v` once per chain of a unit with `n` chains.
_perchain(v::AbstractVector, n) = v
_perchain(v, n) = fill(v, n)

# Concatenate compatible `MarkovChains` (see `_groupkey`) into one. Every table column
# carries over verbatim except the indices into the chains (`hpi`) and `pofs`, an offset
# *within* the latent block, which shift by the running base. A chain's first free point
# always has `lfree == 0`, so a shifted `pofs` never reaches across a chain boundary and
# the coloring scan's carry still resets there.
function MarkovChains(units::AbstractVector{<:MarkovChains})
    u1 = first(units)
    cbase = _bases(map(_nchains, units))
    fbase = _bases(map(_nfree, units))
    free = _cattable(
        [Base.structdiff(u.free, NamedTuple{(:hasfix,)}) for u in units],
        (hpi = cbase, pofs = fbase)
    )
    # A numeric field every chain shares stays a *scalar*, so the arithmetic built from it
    # (σ², 2/τ, the Φ clamp) constant-folds instead of becoming a runtime gather over a
    # constant vector. Only genuinely per-site values pay for the gather.
    numvals = map(keys(u1.numvals)) do k
        vs = reduce(vcat, [_perchain(getproperty(u.numvals, k), _nchains(u)) for u in units])
        return allequal(vs) ? first(vs) : vs
    end
    sites = u1.sites === Val(nothing) ? nothing : reduce((a, b) -> (a..., b...), (_sitetuple(u.sites) for u in units))
    ibase = _bases(map(u -> length(u.inds), units))
    return MarkovChains(
        u1.process, u1.init, Val(sites), u1.param, NamedTuple{keys(u1.numvals)}(numvals),
        reduce(vcat, [u.inds for u in units]), reduce(vcat, [u.ts for u in units]),
        reduce(vcat, map((u, b) -> map(r -> r .+ b, u.spans), units, ibase)),
        reduce(vcat, [u.fixed for u in units]),
        (; free..., hasfix = any(u -> u.free.hasfix, units)),
        _cattable_nested([u.chain for u in units], cbase),
        _cattable_nested([u.fsub for u in units], cbase),
    )
end

_sitetuple(::Val{S}) where {S} = S

_cattable_nested(tabs, cbase) = map(
    k -> _cattable([getproperty(t, k) for t in tabs], (hpi = cbase,)), (open = :open, step = :step)
)

# Partition the chain tuple into maximal runs of compatible chains, preserving order (and
# hence the flat coordinate layout). `IIDChainSpec`s pass through unchanged.
function _build_groups(units::Tuple)
    out = Any[]
    run = MarkovChains[]
    flush!() = (isempty(run) || (push!(out, MarkovChains(run)); run = MarkovChains[]))
    for u in units
        if u isa MarkovChains
            (isempty(run) || isequal(_groupkey(first(run)), _groupkey(u))) || flush!()
            push!(run, u)
        else
            flush!()
            push!(out, u)
        end
    end
    flush!()
    return Tuple(out)
end

# ----- per-point process materialization ---------------------------------------
#
# Each point's own process parameters. Fitted fields of several chains are stacked into a
# length-`nchain` vector (one `concatenate` per field) and gathered per point; fixed fields
# are host constants, which the compiler folds. `setproperties` then rebuilds the process
# exactly as `materialize` does. When every value is a scalar — a single chain, or chains
# sharing every value — the process is built once per walker instead, so it stays
# loop-invariant: its derived quantities hoist out of the loop, and Enzyme accumulates each
# hyperparameter's derivative into a scalar rather than scatter-adding into a vector.
# `_unit_params` is evaluated once per walker and `_unit_process` once per point.

# `vcat(t...)`, NOT `reduce(vcat, t)`. The n-ary form is one `stablehlo.concatenate`; the
# left fold is K-1 pairwise ones, and from the second step on it concatenates a *traced
# array* with a *traced scalar*, which drops off Reactant's fast path into Base's generic
# `typed_vcat` and aborts the trace with "Scalar indexing is disallowed" (Reactant ≥ 0.2.280).
@inline _stackvals(t::Tuple{Any}) = only(t)
@inline _stackvals(t::Tuple) = vcat(t...)

@inline _fieldvals(u::MarkovChains, ::Val{nothing}, hp) = u.numvals
@inline function _fieldvals(u::MarkovChains, ::Val{S}, hp) where {S}
    fitk = _fitkeys(u.process)
    fit = NamedTuple{fitk}(
        map(k -> _stackvals(map(s -> getproperty(getproperty(hp, s), k), S)), fitk)
    )
    return merge(u.numvals, fit)
end

@inline _unit_params(u::MarkovChains, hp) = _shared_process(u.process, _fieldvals(u, u.sites, hp))
@inline _shared_process(p, vals::NamedTuple{F, <:Tuple{Vararg{Number}}}) where {F} = setproperties(p, vals)
@inline _shared_process(p, vals) = vals

@inline _gatherval(v::AbstractVector, j) = rgetindex(v, j)
@inline _gatherval(x::Number, _) = x

@inline _unit_process(::MarkovChains, p::AbstractGaussMarkovProcess, j) = p
@inline _unit_process(u::MarkovChains, vals::NamedTuple, j) =
    setproperties(u.process, map(v -> _gatherval(v, j), vals))
@inline _hpidx(::MarkovChains, free, k) = rgetindex(free.hpi, k)

# ----- exact conditional sampling (OU bridge) -----------------------------------
# rand-only code, never in the logpdf hot path. `_fill_fixed!` zero-initializes, which
# the masked dummy reads of `_free_moments` require.

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

function _rand_chain!(rng::Random.AbstractRNG, x, u::MarkovChains, hp)
    # Wrapped processes have no Gaussian conditionals, so this bridge would silently draw
    # unwrapped values ignoring winding. Each wrapped process registers its own sampler
    # (dispatch, like `MarkovChains{<:WrappedBrownian}` below); this guard makes a missing
    # registration a loud error instead of wrong prior draws.
    is_wrapped(u.process) && throw(
        ArgumentError(
            "no chain sampler registered for the wrapped process $(typeof(u.process).name.name); " *
                "define `_rand_chain!(rng, x, u::MarkovChains{<:$(typeof(u.process).name.name)}, hp)`"
        )
    )
    vals = _unit_params(u, hp)
    free = u.free
    @trace track_numbers = false for k in 1:_nfree(u)
        p = _unit_process(u, vals, _hpidx(u, free, k))
        μ = _process_mean(p)
        m₀, P₀ = initial_moments(u.init, p)
        m, s = _free_moments(p, μ, m₀, P₀, x, free, k)
        rsetindex!(x, m + s * randn(rng, typeof(m)), rgetindex(free.tgt, k))
    end
    return x
end

# Wrapped chains: no Gaussian conditionals, so the generic free-point bridge does not
# apply. Sampling factors into circular random-walk segments around the fixed points:
# a backward walk before the first fixed point (the uniform init makes the time-reversed
# walk a walk again), a forward walk after the last, and, between consecutive fixed
# points, a Brownian bridge whose total 2π winding is drawn exactly from the discrete
# wrapped-normal mixture the wrap induces.
function _sample_winding(rng::Random.AbstractRNG, Δθ, QT)
    # windings beyond ~4 increment sds contribute negligible mass
    W = max(3, ceil(Int, 4 * sqrt(QT) / (2π)))
    # Shift by the smallest exponent (w = 0, since |Δθ| ≤ π) so the dominant weight is
    # exactly 1. Raw exp weights all underflow to 0 once Δθ²/(2QT) ≳ 745 (tight bridges),
    # and the CDF walk over all-zero weights would then deterministically return -W —
    # a spurious multi-turn bridge — instead of the dominant winding.
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

function _rand_chain!(rng::Random.AbstractRNG, x, u::MarkovChains{<:WrappedBrownian}, hp)
    vals = _unit_params(u, hp)
    for j in eachindex(u.spans)
        r = u.spans[j]
        _rand_wrapped_walk!(rng, x, _unit_process(u, vals, j), u.inds[r], u.ts[r], findall(u.fixed[r]))
    end
    return x
end

function _rand_wrapped_walk!(rng::Random.AbstractRNG, x, p, inds, ts, fixpos)
    T = float(eltype(x))
    # A wrapped chain is built only over the stamps where its site is free (`ref_is_gauge`),
    # so a site that is the reference at every one of its stamps has an empty chain. Its
    # values are all reference-fixed and `_fill_fixed!` has already placed them.
    isempty(inds) && return x

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

# Mean-reverting wrapped chains: the same two-sided bridge as the Gaussian processes, but
# with the drift onto the previous chain value taken through the process's own
# `_cond_mean` — the shortest arc for `WrappedOrnsteinUhlenbeck`, the sine drift for
# `VonMisesProcess` — and every drawn value wrapped back onto `(−π, π]`. The innovation is
# drawn normal rather than wrapped normal and the winding of a bridge between two
# reference-fixed points is not enumerated as it is for `WrappedBrownian` above; both are
# the same `σ ≪ π` approximation these processes already document, and this is the
# rand-only path — the log-density is exact.
function _rand_chain!(
        rng::Random.AbstractRNG, x,
        u::MarkovChains{<:Union{WrappedOrnsteinUhlenbeck, VonMisesProcess}}, hp
    )
    vals = _unit_params(u, hp)
    free = u.free
    T = float(eltype(x))
    for k in 1:_nfree(u)
        p = _unit_process(u, vals, _hpidx(u, free, k))
        μ = _process_mean(p)
        m₀, P₀ = initial_moments(u.init, p)
        Φl, Ql = transition_moments(p, free.dtl[k])
        ml = free.mskl[k]
        # prior-from-the-left: the transition from the previous chain value, or the init
        # marginal when this point opens the chain (`ml == 0` masks the drift term off).
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
# A whitened (non-centered) Markov chain's latent coordinates are iid standard variates
# that are *colored* through the triangular map above: each free point is `m + s * z`.
# The map is exact, so it transports the Std spaces (`ascube`/`StdNormal`) with no
# density bookkeeping, while in flat space TV threads the `Σ log s` Jacobian and the
# constrained-space chain logpdf is unchanged. Centered chains keep their raw values as
# coordinates (identity, zero Jacobian) and therefore only support the flat space, where
# the target logpdf is evaluated. IID override chains use their distribution's own node
# per free point (its TV transform in flat space, its exact quantile transport in Std
# spaces).

# Only the angle embedding costs two latent reals per free phase; every other
# parameterization takes one coordinate per free point.
_flat_dim(spec::MarkovChains) =
    spec.param isa AngleEmbedded ? 2 * _nfree(spec) : _nfree(spec)
_flat_dim(spec::IIDChainSpec) = _nfree(spec) * TV.dimension(PT.transport_node(spec.dist, PT.TVFlat()))
_flat_dim(d::GaussMarkovChainDist) = sum(_flat_dim, chaingroups(d); init = 0)

_std_dim(spec::MarkovChains, space) = _nfree(spec)
_std_dim(spec::IIDChainSpec, space) = _nfree(spec) * PT.dimension(PT.transport_node(spec.dist, space))
_std_dim(d::GaussMarkovChainDist, space) = sum(Base.Fix2(_std_dim, space), chaingroups(d); init = 0)

_std_transportable(spec::MarkovChains) = spec.param isa NonCentered && !is_wrapped(spec.process)
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

# --- flat (TransformVariables) chain walkers: thread (logjac, index) over the
# heterogeneous spec tuple with `Base.tail` recursion (TV tuple-transform precedent).

# `index` walks the flat latent vector; `sidx` walks the shared `(av, cv)` staging
# buffers of the batched scan (see `GaussMarkovChainDist.scantgt`). Staging reads `y`
# only at reference-fixed positions (placed by `_fill_fixed!`), never another spec's
# colored values, so all scanning specs stage first and the caller scans+scatters once.
@inline _color_specs_flat!(flag, y, x, index, sidx, av, cv, ::Tuple{}, hp) =
    TV.logjac_zero(flag, eltype(x)), index, sidx
@inline function _color_specs_flat!(flag, y, x, index, sidx, av, cv, specs::Tuple, hp)
    ℓ1, index, sidx = _color_chain_flat!(flag, y, x, index, sidx, av, cv, first(specs), hp)
    ℓ2, index, sidx = _color_specs_flat!(flag, y, x, index, sidx, av, cv, Base.tail(specs), hp)
    return ℓ1 + ℓ2, index, sidx
end

# Shared prologue/epilogue of the flat and Std colorings: the `(a, c)` staging buffers
# for the batched scan over `d.scantgt`, and the single scan+scatter that resolves them.
_scan_buffers(y, d::GaussMarkovChainDist) =
    (similar(y, length(d.scantgt)), similar(y, length(d.scantgt)))
# A non-centered wrapped chain scans the *unwrapped* phase, so its scanned values are
# wrapped onto `(−π, π]` before they are scattered while the real-line chains sharing the
# same buffers are not. The mask makes that one branchless broadcast rather than a second
# scatter over a sub-range (see `GaussMarkovChainDist.scanwrap`).
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

# Wrapped chains: each free phase is embedded as two latent reals through the same
# `AngleTransform` that backs the circular IID priors (`DiagonalVonMises`), so the flat
# coordinates carry no 2π ambiguity. The chain density is still evaluated on the angles by
# `chain_term`, so the exact refant conditioning is untouched, and the transform is
# hyperparameter-independent, so the process is never materialized here.
function _color_chain_wrapped!(flag, y, x, index, spec::MarkovChains)
    free = spec.free
    n = _nfree(spec)
    t = PT.angle_transform()
    # The iterations are independent, but a direct write to the scattered `free.tgt`
    # positions serializes under Reactant's loop raising (a scattered-index scatter never
    # raises; scattered gathers and contiguous writes do). Stage into a contiguous temp
    # and finish with one vectorized `scatter_values!`.
    θs, ℓs = _traced_fill!((similar(y, n), similar(x, n)), (1, 1), n) do k
        θ, ℓk, _ = TV.transform_with(TV.LogJac(), t, x, index + 2 * (k - 1))
        return θ, ℓk
    end
    scatter_values!(y, free.tgt, θs)
    return _logjac_sum(flag, ℓs), index + 2n
end

# Wrapped chains, centered (raw-angle) coordinates: one flat coordinate `φ` per free
# phase, the angle is `φ` wrapped to `(−π, π]`, and the flat density carries the sheet
# weight of `_sheet_logweight` so the lift is *proper* — see the sheet-weight section of
# `processes.jl` for why the reweighting leaves the circular model exactly unchanged.
#
# The sheet anchor of free point `k` is the value the chain density itself conditions on:
# the previous free point's coordinate when the left chain neighbour is free (a
# *contiguous* read of `x`, since free points take one coordinate each in order), the
# reference-fixed value when it is not, and — for the point that opens the chain — the
# window of `_init_window`. All three are sheet-independent, so the sheets of `φₖ` still
# sum to one. Reading the anchor from `x` rather than from the
# colored `y` is what keeps this a pointwise loop: there is no recurrence to resolve, so
# it raises under Reactant exactly like the other stage-then-scatter walkers.
function _color_chain_wrapped_centered!(flag, y, x, index, spec::MarkovChains, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    θs, ℓs = _traced_fill!((similar(y, n), similar(x, n)), (1, 1), n) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        φ = rgetindex(x, index + k - 1)
        ml = rgetindex(free.mskl, k)
        lf = rgetindex(free.lfree, k)
        # Static table rather than a computed index: `free.pofs` points at the previous
        # free point's coordinate, or (masked off) at this point's own.
        prev = rgetindex(x, index + rgetindex(free.pofs, k))
        anchor = lf * prev + (1 - lf) * ml * rgetindex(y, rgetindex(free.lidx, k))
        # Both branches are evaluated and blended by the static mask so the loop body
        # stays single-path; the dummy gap of a chain-opening point is a finite 1 hour,
        # so the unused transition weight is finite rather than NaN.
        wt = _sheet_logweight(p, φ, anchor, rgetindex(free.dtl, k), μ)
        w0 = _init_sheet_logweight(spec.init, p, φ)
        return _wrap_angle(φ), ml * wt + (1 - ml) * w0
    end
    scatter_values!(y, free.tgt, θs)
    return _logjac_sum(flag, ℓs), index + n
end

# Wrapped chains, non-centered coordinates: the coordinate of each free phase is its
# increment onto the previous chain value divided by that step's standard deviation, so the
# unwrapped phase is the running sum `φₖ = φₖ₋₁ + √Qₖ zₖ` and the angle is its wrap. The
# chain's opening point has nothing to take an increment from, so it keeps a raw-angle
# coordinate weighted over `_init_window`, exactly as the centered lift does.
#
# The evaluation is FUSED with `chain_term` (see `_fused_flat`): each stepping point's WN
# transition term there would cancel its `−log WN` sheet weight here exactly, so neither
# is computed — this function contributes `log N(zₖ; 0, 1)` directly (no theta-function
# sums) and `chain_term` only the chain-opening remainder. Only the chain-opening init
# weight keeps a WN factor, one per chain.
#
# `_cond_mean` must be affine in the previous point for the running sum to be a first-order
# affine scan; `_check_param` rejects the processes where it is not.
function _color_chain_wrapped_noncentered!(flag, y, x, index, sidx, av, cv, spec::MarkovChains, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    ℓs = similar(x, n)
    _traced_fill!((av, cv, ℓs), (sidx, sidx, 1), n) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        _, Q = transition_moments(p, rgetindex(free.dtl, k))
        ml = rgetindex(free.mskl, k)
        lf = rgetindex(free.lfree, k)
        z = rgetindex(x, index + k - 1)
        # `a = lfree` resets the carry wherever the previous chain value is not another
        # free point; the fixed value it steps off is then a constant in `c`.
        c = ml * (sqrt(Q) * z + (1 - lf) * rgetindex(y, rgetindex(free.lidx, k))) + (1 - ml) * z
        ℓk = ml * (-(z^2 + log(oftype(float(z), 2π))) / 2) +
            (1 - ml) * _init_sheet_logweight(spec.init, p, z)
        return lf, c, ℓk
    end
    return _logjac_sum(flag, ℓs), index + n, sidx + n
end

# Wrapped chains, non-centered coordinates, for a process whose transition mean carries
# the previous point's sheet but is NOT affine in it (`has_noncentered_lift` without
# `has_affine_lift`, i.e. the sine drift of `VonMisesProcess`). Same coordinates, same
# weights, and the same fused flat target as `_color_chain_wrapped_noncentered!` — each
# stepping point contributes exactly `log N(zₖ; 0, 1)` and a chain opening its raw-angle
# init weight — but the recurrence `θₖ = wrap(m(θₖ₋₁) + √Q zₖ)` is nonlinear, so it cannot
# stage into the batched affine scan and is resolved here as a loop-carried recurrence.
# Equivariance (`m(θ + 2πs) = m(θ) + 2πs`) lets the carry run on the *wrapped* angles: a
# sheet shift of the carry shifts the mean by the same `2π` multiple, which the wrap
# removes, so no unwrapped running phase is needed. Sequential by nature — under Reactant
# a `stablehlo.while` per unit, which is why this parameterization targets scan-level
# chains rather than IntegSeg-scale ones.
function _color_chain_wrapped_recur!(flag, y, x, index, spec::MarkovChains, hp)
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
        # The previous chain value: the carried previous free angle, or the
        # reference-fixed value where the left neighbour is not free (`_fill_fixed!`
        # placed it; at a chain opening `ml` masks the whole step off, and the opening
        # keeps its raw-angle coordinate exactly as in the affine-scan walker).
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

function _color_chain_flat!(flag, y, x, index, sidx, av, cv, spec::MarkovChains, hp)
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
        # A straight copy of a contiguous flat segment to the scattered chain positions:
        # one vectorized scatter (a scattered-write loop would serialize under Reactant).
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
    # The recurrence `y[tgt[k]] = m(y[tgt[k-1]]) + s*x[k]` is affine in the previous
    # free value, so split it: this loop computes the coefficients `(a, c)` — reading
    # only tables, hyperparameters, and fixed values, with contiguous writes, so it
    # raises under Reactant. The recurrence itself is resolved by the caller's single
    # batched `affine_scan!` over all chains.
    ℓs = similar(x, n)
    _traced_fill!((av, cv, ℓs), (sidx, sidx, 1), n) do k
        a, c, s = _unit_affine(spec, vals, y, free, k)
        return a, c + s * rgetindex(x, index + k - 1), log(s)
    end
    return _logjac_sum(flag, ℓs), index + n, sidx + n
end

function _color_chain_flat!(flag, y, x, index, sidx, av, cv, spec::IIDChainSpec, hp)
    dim = TV.dimension(spec.fnode)
    n = _nfree(spec)
    # Independent pointwise transforms; stage contiguously so the loop raises and the
    # scattered write is a single vectorized scatter (see `_color_chain_wrapped!`).
    ys, ℓs = _traced_fill!((similar(y, n), similar(x, n)), (1, 1), n) do k
        yk, ℓk, _ = TV.transform_with(TV.LogJac(), spec.fnode, x, index + (k - 1) * dim)
        return yk, ℓk
    end
    scatter_values!(y, spec.freeinds, ys)
    return _logjac_sum(flag, ℓs), index + n * dim, sidx
end

@inline _whiten_specs_flat!(x, index, y, ::Tuple{}, hp) = index
@inline function _whiten_specs_flat!(x, index, y, specs::Tuple, hp)
    index = _whiten_chain_flat!(x, index, y, first(specs), hp)
    return _whiten_specs_flat!(x, index, y, Base.tail(specs), hp)
end

# Inverse of the centered wrapped lift for a process whose transition mean is itself
# `2π`-periodic (`_periodic_mean`): `_cond_mean` then reads the previous point only through
# its *angle*, so the sheet the coloring's weight favors is pinned by the stored angles
# alone and the unwrap is pointwise. The chain-opening point takes the centre of
# `_init_window`, the same window its weight used.
function _whiten_chain_wrapped_pointwise!(x, index, y, spec::MarkovChains, hp)
    free = spec.free
    vals = _unit_params(spec, hp)
    _traced_fill!(x, index, _nfree(spec)) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        m₀, _ = _init_window(spec.init, p, zero(eltype(x)))
        θ = rgetindex(y, rgetindex(free.tgt, k))
        ml = rgetindex(free.mskl, k)
        Φ, _ = transition_moments(p, rgetindex(free.dtl, k))
        yl = rgetindex(y, rgetindex(free.lidx, k))
        m = ml * _cond_mean(p, yl, Φ, μ) + (1 - ml) * m₀
        return m + _wrap_angle(θ - m)
    end
    return index + _nfree(spec)
end

# Inverse of the centered wrapped lift for a process whose transition mean carries the
# previous point's sheet. The coloring stores wrapped angles, so the coordinates are
# recovered by *unwrapping* the chain: each free point's coordinate is its sheet anchor's
# mean plus the wrapped step onto it, `φₖ = m(φₖ₋₁) + wrap(θₖ − m(φₖ₋₁))`, which places every
# point on the sheet the coloring's weight favors. The mean of the previous *coordinate*
# fixes the sheet, so this is a loop-carried recurrence. `transform ∘ inverse` is the
# identity on wrapped angles; `inverse ∘ transform` is not (the sheet is not recoverable
# from the angle), as for the angle embedding.
function _whiten_chain_wrapped_centered!(x, index, y, spec::MarkovChains, hp)
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
        # sheet anchor of the coloring's weight: the previous coordinate when the left
        # chain point is free, the reference-fixed value when it is not, the init window's
        # centre at a chain opening
        anchor = lf * φp + (1 - lf) * rgetindex(y, rgetindex(free.lidx, k))
        m = ml * _cond_mean(p, anchor, Φ, μ) + (1 - ml) * m₀
        φk = m + _wrap_angle(θ - m)
        rsetindex!(x, φk, index + k - 1)
        φp = φk
    end
    return index + n
end

# Inverse of the non-centered wrapped lift: each coordinate is the wrapped step onto its
# free point's transition mean divided by that step's standard deviation, which places the
# step in `(−π, π]` — the sheet the coloring's weight favors. This is pointwise for
# *every* process the lift admits, the nonlinear sine drift included: the mean is
# `2π`-equivariant, so its value modulo `2π` — all the wrapped step needs — is readable
# from the stored previous angle, with no recurrence to resolve. The chain-opening point's
# coordinate is its angle placed on the sheet nearest its init window's centre, as in the
# coloring. `transform ∘ inverse` is the identity on wrapped angles; `inverse ∘ transform`
# is not (the sheet is not recoverable from the angle), as for the centered lift and the
# angle embedding.
function _whiten_chain_wrapped_noncentered!(x, index, y, spec::MarkovChains, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    _traced_fill!(x, index, n) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        Φ, Q = transition_moments(p, rgetindex(free.dtl, k))
        m₀, _ = _init_window(spec.init, p, zero(eltype(x)))
        θ = rgetindex(y, rgetindex(free.tgt, k))
        ml = rgetindex(free.mskl, k)
        m = _cond_mean(p, rgetindex(y, rgetindex(free.lidx, k)), Φ, μ)
        return ml * (_wrap_angle(θ - m) / sqrt(Q)) + (1 - ml) * (m₀ + _wrap_angle(θ - m₀))
    end
    return index + n
end

function _whiten_chain_flat!(x, index, y, spec::MarkovChains, hp)
    free = spec.free
    if spec.param isa AngleEmbedded
        # Inverse of the angle embedding, `(sin θ, cos θ)` (radius 1); the coloring's
        # `atan` recovers θ exactly, so transform ∘ inverse is the identity on the angles
        # even though inverse ∘ transform is not (the radius is not preserved).
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
            # `_periodic_mean`/`has_affine_lift` are traits on a concrete process type,
            # so these are host-static branches.
            _periodic_mean(spec.process) &&
                return _whiten_chain_wrapped_pointwise!(x, index, y, spec, hp)
            return _whiten_chain_wrapped_centered!(x, index, y, spec, hp)
        end
        _traced_fill!(k -> rgetindex(y, rgetindex(free.tgt, k)), x, index, n)
        return index + n
    end
    if is_wrapped(spec.process)
        return _whiten_chain_wrapped_noncentered!(x, index, y, spec, hp)
    end
    vals = _unit_params(spec, hp)
    _traced_fill!(x, index, _nfree(spec)) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        m₀, P₀ = initial_moments(spec.init, p)
        m, s = _free_moments(p, μ, m₀, P₀, y, free, k)
        return (rgetindex(y, rgetindex(free.tgt, k)) - m) / s
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

# --- Std-space chain walkers: no Jacobian; the latent is converted to a standard
# normal via the space's cdf (exact), mirroring PT's `ScalarTransport`.

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

function _color_chain_std!(y, x, index, sidx, av, cv, spec::MarkovChains, hp, space)
    # centered and wrapped chains are rejected at node construction (`_check_std_transportable`)
    vals = _unit_params(spec, hp)
    free = spec.free
    n = _nfree(spec)
    # Same affine split as the flat coloring (see `_color_chain_flat!`): coefficients
    # staged in raisable loops; the caller's batched `affine_scan!` resolves them.
    _traced_fill!((av, cv), (sidx, sidx), n) do k
        a, c, s = _unit_affine(spec, vals, y, free, k)
        return a, c + s * _to_stdnormal(space, rgetindex(x, index + k - 1))
    end
    return index + n, sidx + n
end

# A latent coordinate of `space` as a standard normal value and back, through the space's cdf.
# In the StdNormal space the latent already is that value, and the map is the identity: the
# cdf/quantile pair would saturate in the tails.
_to_stdnormal(space, z) = PT.space_quantile(PT.StdNormal(), PT._clamp_unit(PT.space_cdf(space, z)))
_to_stdnormal(::PT.StdNormal, z) = z
_from_stdnormal(space, w) = PT.space_quantile(space, PT._clamp_unit(PT.space_cdf(PT.StdNormal(), w)))
_from_stdnormal(::PT.StdNormal, w) = w

function _color_chain_std!(y, x, index, sidx, av, cv, spec::IIDChainSpec, hp, space)
    tin = PT.transport_node(spec.dist, space)
    dim = PT.dimension(tin)
    n = _nfree(spec)
    # Independent pointwise transforms; stage contiguously so the loop raises and the
    # scattered write is a single vectorized scatter (see `_color_chain_wrapped!`).
    ys = _traced_fill!(k -> first(PT.pfwd_step(tin, x, index + (k - 1) * dim)), similar(y, n), 1, n)
    scatter_values!(y, spec.freeinds, ys)
    return index + n * dim, sidx
end

@inline _whiten_specs_std!(x, index, y, ::Tuple{}, hp, space) = index
@inline function _whiten_specs_std!(x, index, y, specs::Tuple, hp, space)
    index = _whiten_chain_std!(x, index, y, first(specs), hp, space)
    return _whiten_specs_std!(x, index, y, Base.tail(specs), hp, space)
end

function _whiten_chain_std!(x, index, y, spec::MarkovChains, hp, space)
    vals = _unit_params(spec, hp)
    free = spec.free
    n = _nfree(spec)
    _traced_fill!(x, index, n) do k
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        m₀, P₀ = initial_moments(spec.init, p)
        m, s = _free_moments(p, μ, m₀, P₀, y, free, k)
        return _from_stdnormal(space, (rgetindex(y, rgetindex(free.tgt, k)) - m) / s)
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

# Flat node for the fixed-hyperparameter case; wrapped in `InstrumentTransform` by
# `ObservedArrayPrior`'s generic node, so it returns the plain full-length vector.
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

# One sample of a hierarchical instrument prior: the gains plus the fitted process
# hyperparameters. Every method that recognizes the shape by dispatch (`siteparams`,
# `caltable`, the Makie `prepare_caltable`) uses this alias, so the shape is written
# down once.
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

#    ObservedHierarchicalArrayPrior
#
#The observed prior produced by an `ArrayPrior` whose [`GaussMarkovSitePrior`](@ref)
#processes have fitted hyperparameters. Samples are NamedTuples
#`(params = SiteArray, hyperparams = NamedTuple)` and the log-density is the exact chain
#density of the parameters given the hyperparameters plus the hyperprior. Internal type.
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

# Hierarchical transport: the hyperparameter coordinates come *first* in the latent
# vector so the node can materialize the processes before coloring the chains. The
# node's shape depends only on the sitemap, the fixed indices, and *which* fields are
# fitted — never on hyperparameter values — preserving the static-transport invariant.
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

# Chain times in *absolute* UTC hours — the same clock as the data's `Ti` field — with
# multi-day tracks continued past 24 h through the `mjd` field. The origin is the whole
# array's first day, shared by every chain, so a chain's times are directly comparable
# both to the data and between sites; times relative to each chain's own first point
# would make a chain's `ts` silently unmatchable against `Ti` and incomparable across
# sites. Only *gaps* reach the process — `transition_moments`, the bridge, and the
# `Δ0` of `_first_point_term` are all differences — so the origin never enters the
# density; this is purely an interface property.
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

# Split a chain's time stamps into the ones it is built over and, among those, the ones it
# is conditioned on. A process whose reference values are a gauge (`ref_is_gauge`) drops
# them from the chain entirely, so the free points step off each other over the composed
# gap. A `FixedInit` chain's opening point is exempt: that value states where the chain
# starts, which is a statement about the process and not a gauge, so it is retained and
# stays conditioned on. Everything downstream reads the chain through `inds`/`ts` and the
# tables built from them, so this is the only place the distinction appears.
function _gauge_split(process, init, n::Int, fixedpos)
    ref_is_gauge(process) || return collect(1:n), collect(fixedpos)
    opens = init isa FixedInit && !isempty(fixedpos) && first(fixedpos) == 1
    drop = Set(opens ? fixedpos[2:end] : fixedpos)
    return filter(!in(drop), 1:n), opens ? [1] : Int[]
end

function _chainspec(sp::GaussMarkovSitePrior, site, inds, ts, fixedpos)
    hpk = isempty(hyperprior(sp.process)) ? nothing : site
    keep, kept = _gauge_split(sp.process, sp.init, length(inds), fixedpos)
    return MarkovChains(sp.process, sp.init, Val(hpk), inds[keep], ts[keep], kept, sp.param)
end

function _chainspec(sp::IIDSitePrior, site, inds, ts, fixedpos)
    free = isempty(fixedpos) ? copy(inds) : inds[setdiff(eachindex(inds), fixedpos)]
    return IIDChainSpec(sp.dist, free)
end

# Working float precision of the chains, derived from the process (and IID-override)
# parameters so that the sample element type follows the inputs rather than a hardcoded
# `Float64`. Fitted (Distribution) fields contribute their `eltype`.
_site_paramtype(sp::GaussMarkovSitePrior) = promote_type(_paramtype(sp.process), _init_paramtype(sp.init))
_site_paramtype(sp::IIDSitePrior) = float(eltype(sp.dist))
function _working_type(site_dists::NamedTuple)
    T = mapreduce(_site_paramtype, promote_type, values(site_dists); init = Union{})
    return T === Union{} ? Float64 : T
end

function build_markov_observed(d::ArrayPrior, site_dists::NamedTuple, smap::SiteLookup, t)
    d.centroid_station === nothing ||
        throw(ArgumentError("centroid_station is not supported with GaussMarkovSitePrior"))
    # The chain-machinery path has no first-stamp pin for IID sites; refusing here beats
    # silently dropping the pin the prior asked for.
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
    finds, vals = reference_indices(t, smap, d.refant)
    fixedinds = collect(Int, finds)
    fixedvals = vals === nothing ? T[] : collect(T, vals)

    # FixedInit chains are exactly conditioned to the init value at each chain's first
    # time, on top of any reference-fixed indices (see the `init` docs of
    # `GaussMarkovSitePrior`). A referencing scheme that already fixes a chain's first
    # point must agree with the init value; a silent precedence would break the
    # documented FixedInit invariant.
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

    # Set membership: `in(fixedinds)` over a Vector is O(|fixedinds|) per element, which
    # at IntegSeg scale makes this loop quadratic in the number of time stamps.
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

    # Assemble the merged hyperprior: every site gets its own hyperparameters, nested
    # under the site name. The process in the ArrayPrior default is a per-site template,
    # so its fitted fields are replicated for each site that uses it.
    hp = (;)
    for (s, sp) in pairs(site_dists)
        sp isa GaussMarkovSitePrior || continue
        h = hyperprior(sp.process)
        isempty(h) && continue
        hp = merge(hp, NamedTuple{(s,)}((h,)))
    end

    dists = GaussMarkovChainDist(chains, fixedinds, fixedvals, length(smap.times))
    isempty(hp) && return ObservedArrayPrior(dists, smap)
    return ObservedHierarchicalArrayPrior(dists, NamedDist(hp), smap)
end
