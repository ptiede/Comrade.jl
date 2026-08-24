export GaussMarkovSitePrior, siteparams

"""
    GaussMarkovSitePrior(seg::TimeSegmentation, process::AbstractGaussMarkovProcess; init = StationaryInit(), centered = false)

A site prior that is correlated in time, following the continuous-time Gauss-Markov
`process` (e.g. [`OrnsteinUhlenbeck`](@ref)) sampled at the times implied by the
segmentation `seg`. The log-density factorizes over each site's time chain using the
process's exact discretization, so the cost is `O(n)` in the number of times and the
irregular time sampling of a VLBI observation is handled exactly. When the observation
has multiple frequency channels each (site, channel) pair is an independent chain.

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

The default value of `centered` depends on the process, because the trade-off does.

For a **real-line process** (`OrnsteinUhlenbeck`, `BrownianMotion`) the default is
`centered = false`, i.e. the prior is *whitened* (non-centered): each free point's latent
coordinate is an iid standard variate that is colored through the chain's exact
conditional moments. This removes the hyperparameter-gain funnel when hyperparameters are
fitted and makes the transport to the Std spaces exact, so both `asflat` and `ascube`
work. Set `centered = true` to instead use the gain values themselves as coordinates,
which can mix better when every point is strongly data-constrained; centered priors
support `asflat` only.

For a **wrapped (circular) process** ([`WrappedBrownian`](@ref),
[`WrappedOrnsteinUhlenbeck`](@ref)) there is no whitening —
a wrapped chain has no Gaussian conditionals — so the choice is between raw angles and an
embedding, and the default is `centered = true`: one coordinate per free phase, the angle
itself. A wrapped chain's density is exactly `2π`-periodic, so simply *calling* a real
coordinate the angle would make the flat target improper — an infinite lattice of
identical modes, `2π` apart in every coordinate, which a sampler drifts through and which
wrecks warmup variance adaptation. The raw-angle lift therefore carries a *sheet weight*:
each free phase is reweighted by `N(δ; 0, Q) / WN(δ; Q)`, where `δ` is its step onto the
previous chain value and `Q` that step's variance. Those weights sum to exactly one over
the `2π` sheets, so the circular model is untouched — no mass is added, moved, or lost —
while the flat target becomes the proper *unwrapped* Gaussian random walk, in which a
`2π` step costs `exp(−(2π)²/2Q)`.

`centered = false` on a wrapped process selects the alternative: each free phase is
embedded as two latent reals through the same angle transform used by `DiagonalVonMises`.
The embedding is the one parameterization with no sheet structure at all (a phase wrap is
a continuous winding of the latent point), so it is the safe fallback if a fit ever does
show sheet-related pathology — at the cost of twice the phase dimension and a per-pair
rotated-ellipse geometry a diagonal mass matrix cannot capture. Wrapped chains support
`asflat` only in either form.

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
wrapped process removes. In all cases the `phase=true` cumulative reparameterization of
`ArrayPrior` is superseded by this prior: `phase` must be `false`. Reference stations
(`refant`) are handled by exact conditioning of the chain on the fixed values, which
works with scattered fixed indices such as those produced by `SEFDReference`.

The initial distribution of each chain — how its first time stamp is treated — is set by
`init` (see [`AbstractInitialPrior`](@ref)). [`StationaryInit`](@ref) (the default)
starts the chain in the process's stationary marginal. [`GaussianInit`](@ref) uses an
explicit `N(μ0, σ0²)` first marginal, e.g. an intentionally wide one so the starting
value is essentially fit freely. [`FixedInit`](@ref) conditions the chain to an exact
value at its first time, the continuous-time analogue of the cumulative `phase=true`
parameterization: when the chain models an unwrapped phase fluctuation on top of a
separate (circular) offset term, `FixedInit(0.0)` removes the level freedom that would
otherwise trade off against the offset and, through the likelihood's `2π` periodicity,
produce spurious posterior modes at `2π`-shifted levels. A referencing scheme that
already fixes a chain's first point must agree with the `FixedInit` value; a conflict is
an error.

# Example

```julia
lg ~ ArrayPrior(GaussMarkovSitePrior(IntegSeg(), OrnsteinUhlenbeck(σ = VLBIExponential(0.1), τ = VLBIExponential(2.0))))
gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 3.0, τ = 1.0)); refant = SEFDReference(0.0))
## phase-fluctuation term meant to be combined with a circular per-track offset
dgp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = 1.0); init = FixedInit(0.0)); refant = SEFDReference(0.0))
```

!!! note
    When any hyperparameters are fitted the example above uses the `VLBI*` distributions
    (`VLBIExponential`, `VLBIInverseGamma`, ...). These accept a scalar traced argument and
    so work on both the CPU and GPU/`Reactant` backends; the plain `Distributions.jl` types
    only work on the CPU.
"""
struct GaussMarkovSitePrior{S <: TimeSegmentation, P <: AbstractGaussMarkovProcess, I <: AbstractInitialPrior} <: AbstractSitePrior
    seg::S
    process::P
    init::I
    centered::Bool
end

# Wrapped (circular) chains default to the centered raw-angle coordinates; real-line
# chains default to the whitened ones. The two defaults differ because the trade-off
# does: whitening a real-line chain removes the hyperparameter funnel and buys the Std
# transports, while a wrapped chain has no whitening (no Gaussian conditionals) and its
# alternative — the angle embedding — costs twice the dimension for a geometry a diagonal
# mass matrix cannot fit.
function GaussMarkovSitePrior(
        seg::TimeSegmentation, process::AbstractGaussMarkovProcess;
        init::AbstractInitialPrior = StationaryInit(), centered = is_wrapped(process)
    )
    _check_init(init, process)
    return GaussMarkovSitePrior(seg, process, init, centered)
end

needs_chain_machinery(::GaussMarkovSitePrior) = true

# ----- per-(site, freq) chain specifications --------------------------------
#
# Built once at `set_array` time; everything in a spec is constant during sampling
# (indices, times, and the process *template* whose fitted fields are placeholders).

struct MarkovChainSpec{P <: AbstractGaussMarkovProcess, I <: AbstractInitialPrior, K, V <: AbstractVector{<:Integer}, T <: AbstractVector, F <: NamedTuple, G <: NamedTuple}
    process::P     # template; fitted fields are Distribution placeholders
    init::I        # initial distribution of the chain's first time stamp
    hpsel::Val{K}  # the site whose hyperparameters this chain reads, or nothing if fully fixed
    inds::V        # flat indices into the full parameter vector, ascending in time
    ts::T          # chain times in hours (strictly increasing)
    centered::Bool # coordinates are raw values (centered) or whitened standard variates
    free::F        # static per-free-point tables (see `_free_point_tables`)
    fsub::G        # inds/ts restricted to the fixed points, materialized (views do not trace)
end

function MarkovChainSpec(
        process::AbstractGaussMarkovProcess, init::AbstractInitialPrior, hpsel::Val, inds, ts, fixedpos, centered::Bool
    )
    return MarkovChainSpec(
        process, init, hpsel, inds, ts, centered,
        _free_point_tables(inds, ts, fixedpos),
        (inds = inds[fixedpos], ts = ts[fixedpos]),
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

EnzymeRules.inactive_type(::Type{<:MarkovChainSpec}) = true
EnzymeRules.inactive_type(::Type{<:IIDChainSpec}) = true

struct MarkovGroup{P <: AbstractGaussMarkovProcess, I <: AbstractInitialPrior, S, N <: NamedTuple, F <: NamedTuple, C <: NamedTuple}
    process::P     # template of the *first* member; fitted fields are Distribution placeholders
    init::I
    sites::Val{S}  # member sites in group order, or `nothing` when no field is fitted
    centered::Bool
    numvals::N     # per-chain host vectors of the numeric (non-fitted) process fields
    free::F        # concatenated `_free_point_tables`, plus `hpi` (point -> member index)
    chain::C       # concatenated per-time-stamp tables for the chain log-density
    fsub::C        # `chain`, restricted to the reference-fixed time stamps
end

EnzymeRules.inactive_type(::Type{<:MarkovGroup}) = true

# A single chain and a batched group expose the same surface to every walker: `process`,
# `init`, `centered`, `free` and `_nfree`. See the batched-group section below.
const ChainUnit = Union{MarkovChainSpec, MarkovGroup}

@inline _selecthp(hp::NamedTuple, ::Val{nothing}) = (;)
@inline _selecthp(hp::NamedTuple, ::Val{K}) where {K} = getproperty(hp, K)

# ----- the chain log-density -------------------------------------------------

# Log-density of the chain's first point under the init marginal propagated a gap Δ0
# from the chain start (Δ0 = 0 for the full chain; the first fixed time for the
# conditioning subtraction in `chain_term`). `z` is the `iszero(Δ0)` indicator that keeps
# `marginal_moments` branchless: it is a plain host comparison for a single chain and a
# static table read for a batched group, so one family serves both. For FixedInit the first
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

# Exact conditioning on the reference-fixed values: the restriction of an order-1 Markov
# process to a subset of times is Markov with the composed transition law over the larger
# gaps, and its first marginal is the init marginal propagated to the subset's first time
# (`marginal_moments`). Hence log p(free | fixed) = ℓ(all) − ℓ(fixed subset). The
# subtracted term depends on the hyperparameters even though the fixed values are
# constants — dropping it would bias the hyperparameter posterior.
@inline function chain_term(spec::MarkovChainSpec, x, hp)
    # A fully reference-fixed chain conditions on itself: the two logpdf passes below are
    # identical and cancel exactly for any hyperparameters, so skip both (host-static
    # branch — the index vectors are concrete).
    length(spec.fsub.inds) == length(spec.inds) && return zero(eltype(x))
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    ℓ = _gm_chain_logpdf(p, spec.init, x, spec.inds, spec.ts)
    isempty(spec.fsub.inds) && return ℓ
    return ℓ - _gm_chain_logpdf(p, spec.init, x, spec.fsub.inds, spec.fsub.ts, first(spec.fsub.ts) - first(spec.ts))
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
#[`GaussMarkovSitePrior`](@ref): a set of independent per-(site, frequency) Gauss-Markov
#chains (plus IID terms for `IIDSitePrior` overrides), exactly conditioned on the
#reference-fixed indices. This is an internal type constructed by `ObservedArrayPrior`.
struct GaussMarkovChainDist{C <: NamedTuple, G <: Tuple, I <: AbstractVector{<:Integer}, F} <: Dists.ContinuousMultivariateDistribution
    chains::C
    # `chains` batched into `MarkovGroup`s (plus any `IIDChainSpec`s, unchanged), in spec
    # order, used for the dimension/trait queries and by the traced walkers (see
    # `_walk_units`). `chains` itself is what the host walkers and the sampler iterate.
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
end

_scan_tgt(spec::ChainUnit) =
    (is_wrapped(spec.process) || spec.centered) ? Int[] : spec.free.tgt
_scan_tgt(::IIDChainSpec) = Int[]

function GaussMarkovChainDist(chains::NamedTuple, fixedinds, fixedvals, len::Int)
    groups = _build_groups(values(chains))
    scantgt = reduce(vcat, map(_scan_tgt, collect(groups)); init = Int[])
    return GaussMarkovChainDist(chains, groups, fixedinds, fixedvals, len, scantgt)
end

# The batched units (see `_build_groups`).
chaingroups(d::GaussMarkovChainDist) = d.groups
EnzymeRules.inactive(::typeof(chaingroups), args...) = nothing

# Which units the per-evaluation walkers iterate. The batched group is the concatenation of
# its chains *in order*, so both choices produce the same flat layout, the same `scantgt`,
# and the same numbers — they differ only in cost model, and the two cost models disagree:
#
#   - On an accelerator every chain costs a kernel launch, and the arithmetic is free at
#     these sizes. Batching makes the launch count O(1) in the number of chains instead of
#     O(nchains) (measured post-XLA: 670 -> 124 kernels at 7 sites, and flat rather than
#     linear as sites are added).
#   - On the host the loops are plain Julia and launches do not exist. The per-chain form
#     keeps the process loop-invariant, so `σ²`/`stationary_moments` hoist out of the loop
#     and Enzyme accumulates `∂/∂σ` into a *scalar* shadow. The group reads each point's
#     parameters through a gather, which makes that a scatter-add into a length-nchain
#     shadow and costs Enzyme a per-iteration allocation: measured 3.4x slower gradients
#     (34 -> 118 μs, 21 -> 76 allocs) for a fitted-hyperparameter OU chain.
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

_nfree(spec::MarkovChainSpec) = length(spec.free.tgt)
_nfree(spec::IIDChainSpec) = length(spec.freeinds)

# ----- batched chain groups ----------------------------------------------------
#
# Every walker below used to run once per (site, frequency) chain: its own staged loop,
# its own `scatter_values!`, its own reduction and — with fitted hyperparameters — its own
# scalar pipeline through `materialize`. On an accelerator each of those is a separate
# kernel over a few dozen elements, so the cost of the prior scaled with the *number of
# chains* rather than with the amount of arithmetic. (Measured on the post-XLA gradient
# module at fixed total point count: ~4 extra kernels per Gaussian chain and ~37 per
# wrapped chain, with the arithmetic itself unchanged.)
#
# A `MarkovGroup` is the concatenation of a *run of adjacent compatible chains*: the same
# static tables of `_free_point_tables`, laid out chain after chain, plus a per-point index
# into the group's per-chain process parameters. The loop bodies are the per-chain ones
# unchanged — they simply run once over every point of every chain, so Reactant raises
# each into a single set of vectorized ops instead of one set per chain. Only *adjacent*
# specs are grouped, so the flat coordinate layout is bit-for-bit the one the per-chain
# walkers produced.
#
# Chains may be grouped when they agree on everything the loop bodies dispatch on: the
# process *type* (which also fixes which fields are fitted), the initial prior (by value —
# `FixedInit`/`GaussianInit` carry numbers the moments depend on), and `centered`. They
# need *not* agree on the numeric field values of their process templates: every field is
# read per point from a length-`nchain` vector, so per-site overrides that only change a
# fixed σ or τ still share one group. `hasfix` need not agree either — a chain with no
# reference-fixed points carries `mskf = 0`, and `_bridge_weights(Q₁, 0, 1)` collapses
# exactly to the one-sided conditional, so mixing costs a masked-out bridge term and
# changes no value.


_nfree(g::MarkovGroup) = length(g.free.tgt)

# Chains are grouped iff they agree on what the loop bodies dispatch on. `init` and
# `centered` compare by value; the process compares by *type*, which also fixes which
# fields are `Distribution`s and hence whether `hpsel` is a site or `nothing`.
_groupkey(spec::MarkovChainSpec) = (typeof(spec.process), spec.init, spec.centered)

# Field names of `p` that are fitted (Distribution) and that are fixed (numeric).
_fitkeys(p::AbstractGaussMarkovProcess) = keys(hyperprior(p))
_numkeys(p::AbstractGaussMarkovProcess) = filter(!Base.Fix2(in, _fitkeys(p)), propertynames(p))

# Concatenate the members' `_free_point_tables`. Every column carries over verbatim except
# `pofs`, which is an offset *within* a chain's latent block and so shifts by the running
# base — the rule itself stays in `_free_point_tables`. `hpi` (point -> member) and the
# `hasfix` reduction are the only genuinely new columns. A chain's first free point always
# has `lfree == 0`, so a shifted `pofs` never reaches across a chain boundary and the
# coloring scan's carry still resets there.
function _group_free(specs)
    ks = (:tgt, :lidx, :dtl, :mskl, :lfree, :fidx, :dtf, :mskf)
    cols = NamedTuple{ks}(map(k -> reduce(vcat, (getproperty(sp.free, k) for sp in specs)), ks))
    ns = [length(sp.free.tgt) for sp in specs]
    bases = cumsum(ns) .- ns
    pofs = reduce(vcat, (sp.free.pofs .+ b for (sp, b) in zip(specs, bases)))
    hpi = reduce(vcat, (fill(j, n) for (j, n) in enumerate(ns)))
    return (; cols..., pofs, hpi, hasfix = any(sp -> sp.free.hasfix, specs))
end

# Per-time-stamp tables for the batched chain log-density. Each entry carries the flat
# index of its point and of the point before it *in the same chain*, the gap between them,
# and a mask marking the chain-opening entries, whose transition term is replaced by the
# initial-marginal term. A chain-opening entry reads itself as its predecessor over a
# dummy 1 hour gap, so the masked-off transition stays finite rather than `0 * NaN`.
# `Δ0` propagates the initial marginal to the first time stamp of a subset (see
# `chain_term`); `z0` is its host-computed `iszero` indicator, which keeps the
# `marginal_moments` shortcut branchless now that `Δ0` is a per-point table read.
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
    # A numeric field whose value every member shares stays a *scalar*, so the arithmetic
    # built from it (σ², 2/τ, the Φ clamp) constant-folds exactly as it did per chain
    # instead of becoming a runtime gather over a constant vector. Only genuinely per-site
    # overrides pay for the gather.
    numvals = NamedTuple{_numkeys(s1.process)}(
        map(_numkeys(s1.process)) do k
            vs = [getproperty(sp.process, k) for sp in specs]
            return allequal(vs) ? first(vs) : vs
        end
    )
    # A fully reference-fixed chain conditions on itself: its two log-density passes are
    # identical for any hyperparameters and cancel, so drop it from both tables rather
    # than emitting the cancelling work (the per-chain `chain_term` short-circuits it).
    live = [(j, sp) for (j, sp) in enumerate(specs) if length(sp.fsub.inds) != length(sp.inds)]
    chain = _group_chain(((j, (sp.inds, sp.ts, zero(T))) for (j, sp) in live), T)
    # Only chains that actually have reference-fixed points contribute the conditioning
    # subtraction; the rest have nothing to condition on.
    fsub = _group_chain(
        (
            (j, (sp.fsub.inds, sp.fsub.ts, first(sp.fsub.ts) - first(sp.ts)))
                for (j, sp) in live if !isempty(sp.fsub.inds)
        ),
        T
    )
    return MarkovGroup(
        s1.process, s1.init, Val(sites), s1.centered, numvals,
        _group_free(specs), chain, fsub
    )
end

_hpsite(spec::MarkovChainSpec{<:Any, <:Any, K}) where {K} = K

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
# The one thing a batched loop needs that a per-chain loop does not: each point's own
# process parameters. Fitted fields are stacked into a length-`nchain` vector (one
# `concatenate` per field) and gathered per point; fixed fields are host constants, which
# the compiler folds. `setproperties` then rebuilds the process exactly as `materialize`
# does, so every downstream formula is byte-identical to the per-chain path.

# `vcat(t...)`, NOT `reduce(vcat, t)`. The n-ary form is one `stablehlo.concatenate`; the
# left fold is K-1 pairwise ones, and from the second step on it concatenates a *traced
# array* with a *traced scalar*, which drops off Reactant's fast path into Base's generic
# `typed_vcat` and aborts the trace with "Scalar indexing is disallowed" (Reactant ≥ 0.2.280).
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

# A field the whole group shares is stored as a scalar and read straight through (see the
# `numvals` comment); only a per-site one is gathered.
@inline _gatherval(v::AbstractVector, j) = rgetindex(v, j)
@inline _gatherval(x::Number, _) = x

@inline function _point_process(g::MarkovGroup, vals::NamedTuple{F}, j) where {F}
    return setproperties(g.process, NamedTuple{F}(map(v -> _gatherval(v, j), values(vals))))
end

# A single chain and a batched group expose the same surface to every walker below —
# `process`, `init`, `centered`, `free`, `_nfree` — so the walkers have one body each. The
# only difference is where a point's process parameters come from: a single chain
# materializes once up front, a group gathers per point. `_unit_params` is evaluated once
# per walker, `_unit_process` once per point, and both collapse to the old code for a
# single chain.

@inline _unit_params(spec::MarkovChainSpec, hp) = materialize(spec.process, _selecthp(hp, spec.hpsel))
@inline _unit_process(::MarkovChainSpec, p, j) = p
@inline _hpidx(::MarkovChainSpec, free, k) = 1

@inline _unit_params(g::MarkovGroup, hp) = _group_fieldvals(g, hp)
@inline _unit_process(g::MarkovGroup, vals, j) = _point_process(g, vals, j)
@inline _hpidx(::MarkovGroup, free, k) = rgetindex(free.hpi, k)

# ----- the batched chain log-density --------------------------------------------

function _group_chain_logpdf(g::MarkovGroup, tab, x, vals)
    ℓ = zero(eltype(x))
    @trace track_numbers = false for i in 1:length(tab.inds)
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
    vals = _group_fieldvals(g, hp)
    ℓ = _group_chain_logpdf(g, g.chain, x, vals)
    isempty(g.fsub.inds) && return ℓ
    return ℓ - _group_chain_logpdf(g, g.fsub, x, vals)
end

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

function _rand_chain!(rng::Random.AbstractRNG, x, spec::MarkovChainSpec, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    # Wrapped processes have no Gaussian conditionals, so this bridge would silently draw
    # unwrapped values ignoring winding. Each wrapped process registers its own sampler
    # (dispatch, like `MarkovChainSpec{<:WrappedBrownian}` below); this guard makes a
    # missing registration a loud error instead of wrong prior draws.
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

function _rand_chain!(rng::Random.AbstractRNG, x, spec::MarkovChainSpec{<:WrappedBrownian}, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    inds = spec.inds
    ts = spec.ts
    T = float(eltype(x))

    # positions of the fixed points, in chain order (`fsub.inds` preserves it)
    fx = spec.fsub.inds
    fixpos = Vector{Int}(undef, length(fx))
    j = 1
    @inbounds for i in eachindex(inds)
        j > length(fx) && break
        if inds[i] == fx[j]
            fixpos[j] = i
            j += 1
        end
    end

    step_walk(prev, Δt) = _wrap_angle(prev + sqrt(_increment_var(p, Δt)) * randn(rng, T))

    if isempty(fixpos)
        x[inds[1]] = T(π) * (2rand(rng, T) - 1)
        @inbounds for i in 2:length(inds)
            x[inds[i]] = step_walk(x[inds[i - 1]], ts[i] - ts[i - 1])
        end
        return x
    end

    # backward walk before the first fixed point
    @inbounds for i in (fixpos[1] - 1):-1:1
        x[inds[i]] = step_walk(x[inds[i + 1]], ts[i + 1] - ts[i])
    end
    # bridges between consecutive fixed points, with the winding drawn exactly
    @inbounds for s in 1:(length(fixpos) - 1)
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
    @inbounds for i in (fixpos[end] + 1):length(inds)
        x[inds[i]] = step_walk(x[inds[i - 1]], ts[i] - ts[i - 1])
    end
    return x
end

# Mean-reverting wrapped chains: the same two-sided bridge as the Gaussian processes, but
# with every deviation from the circular mean taken along the shortest arc (matching
# `_cond_mean`) and every drawn value wrapped back onto `(−π, π]`. The innovation is drawn
# normal rather than wrapped normal and the winding of a bridge between two reference-fixed
# points is not enumerated as it is for `WrappedBrownian` above; both are the same `σ ≪ π`
# approximation `WrappedOrnsteinUhlenbeck` already documents, and this is the rand-only
# path — the log-density is exact.
function _rand_chain!(rng::Random.AbstractRNG, x, spec::MarkovChainSpec{<:WrappedOrnsteinUhlenbeck}, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ = _process_mean(p)
    m₀, P₀ = initial_moments(spec.init, p)
    free = spec.free
    T = float(eltype(x))
    @inbounds for k in 1:_nfree(spec)
        Φl, Ql = transition_moments(p, free.dtl[k])
        ml = free.mskl[k]
        # prior-from-the-left: the transition from the previous chain value, or the init
        # marginal when this point opens the chain (`ml == 0` zeroes the Φ₁ term).
        Φ₁ = ml * Φl
        Q₁ = ml * Ql + (1 - ml) * P₀
        b₁ = Φ₁ * _wrap_angle(x[free.lidx[k]] - μ) + (1 - ml) * (m₀ - μ)
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

# Non-centered wrapped chains embed each free phase as two latent reals (AngleTransform);
# centered ones use the raw angles.
_flat_dim(spec::ChainUnit) =
    (is_wrapped(spec.process) && !spec.centered) ? 2 * _nfree(spec) : _nfree(spec)
_flat_dim(spec::IIDChainSpec) = _nfree(spec) * TV.dimension(PT.transport_node(spec.dist, PT.TVFlat()))
_flat_dim(d::GaussMarkovChainDist) = sum(_flat_dim, chaingroups(d); init = 0)

_std_dim(spec::ChainUnit, space) = _nfree(spec)
_std_dim(spec::IIDChainSpec, space) = _nfree(spec) * PT.dimension(PT.transport_node(spec.dist, space))
_std_dim(d::GaussMarkovChainDist, space) = sum(Base.Fix2(_std_dim, space), chaingroups(d); init = 0)

_std_transportable(spec::ChainUnit) = !spec.centered && !is_wrapped(spec.process)
_std_transportable(::IIDChainSpec) = true
function _check_std_transportable(d::GaussMarkovChainDist)
    all(_std_transportable, chaingroups(d)) || throw(
        ArgumentError(
            "This GaussMarkovSitePrior cannot be transported to the Std spaces " *
                "(ascube/StdNormal): the centered parameterization needs the target " *
                "log-density, which those transports never evaluate, and a wrapped " *
                "(circular) chain has no measure-preserving map to the line. Use asflat() " *
                "(or, for centered chains, the default whitened form)."
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
function _resolve_scan!(y, d::GaussMarkovChainDist, av, cv)
    isempty(d.scantgt) || scatter_values!(y, d.scantgt, affine_scan!(av, cv))
    return y
end

function _color_all_flat!(flag, y, x, index, d::GaussMarkovChainDist, hp)
    av, cv = _scan_buffers(y, d)
    ℓ, index, _ = _color_specs_flat!(flag, y, x, index, 1, av, cv, _walk_units(d, x), hp)
    _resolve_scan!(y, d, av, cv)
    return ℓ, index
end

# Skip the `log` on the value-only path; dispatch is on the compile-time flag, so the
# traced loop body stays single-path.
@inline _maybe_logjac(::TV.NoLogJac, s) = zero(s)
@inline _maybe_logjac(::TV.LogJac, s) = log(s)

# Under NoLogJac the per-point logjac from `TV.transform_with` is the `NoLogJac` sentinel,
# not a number; keep the traced accumulator real on both paths (cf. `_maybe_logjac`).
@inline _acc_logjac(::TV.NoLogJac, ℓk, θ) = zero(θ)
@inline _acc_logjac(::TV.LogJac, ℓk, θ) = ℓk

# Wrapped chains: each free phase is embedded as two latent reals through the same
# `AngleTransform` that backs the circular IID priors (`DiagonalVonMises`), so the flat
# coordinates carry no 2π ambiguity. The chain density is still evaluated on the angles by
# `chain_term`, so the exact refant conditioning is untouched, and the transform is
# hyperparameter-independent, so the process is never materialized here.
function _color_chain_wrapped!(flag, y, x, index, spec::ChainUnit)
    free = spec.free
    n = _nfree(spec)
    t = PT.angle_transform()
    ℓ = zero(eltype(x))
    # The iterations are independent, but a direct write to the scattered `free.tgt`
    # positions serializes under Reactant's loop raising (a scattered-index scatter never
    # raises; scattered gathers and contiguous writes do). Stage into a contiguous temp
    # and finish with one vectorized `scatter_values!`.
    θs = similar(y, n)
    @trace track_numbers = false for k in 1:n
        θ, ℓk, _ = TV.transform_with(flag, t, x, index + 2 * (k - 1))
        rsetindex!(θs, θ, k)
        ℓ += _acc_logjac(flag, ℓk, θ)
    end
    scatter_values!(y, free.tgt, θs)
    flag isa TV.NoLogJac && return TV.logjac_zero(flag, eltype(x)), index + 2n
    return ℓ, index + 2n
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
function _color_chain_wrapped_centered!(flag, y, x, index, spec::ChainUnit, hp)
    free = spec.free
    n = _nfree(spec)
    vals = _unit_params(spec, hp)
    ℓ = zero(eltype(x))
    θs = similar(y, n)
    @trace track_numbers = false for k in 1:n
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
        ℓ += ml * wt + (1 - ml) * w0
        rsetindex!(θs, _wrap_angle(φ), k)
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
    if spec.centered
        if is_wrapped(spec.process)
            ℓc, index = _color_chain_wrapped_centered!(flag, y, x, index, spec, hp)
            return ℓc, index, sidx
        end
        # A straight copy of a contiguous flat segment to the scattered chain positions:
        # one vectorized scatter (a scattered-write loop would serialize under Reactant).
        scatter_values!(y, free.tgt, view(x, index:(index + n - 1)))
        return ℓ0, index + n, sidx
    end
    if is_wrapped(spec.process)
        ℓw, index = _color_chain_wrapped!(flag, y, x, index, spec)
        return ℓw, index, sidx
    end
    vals = _unit_params(spec, hp)
    ℓ = zero(eltype(x))
    # The recurrence `y[tgt[k]] = m(y[tgt[k-1]]) + s*x[k]` is affine in the previous
    # free value, so split it: this loop computes the coefficients `(a, c)` — reading
    # only tables, hyperparameters, and fixed values, with contiguous writes, so it
    # raises under Reactant. The recurrence itself is resolved by the caller's single
    # batched `affine_scan!` over all chains.
    @trace track_numbers = false for k in 1:n
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        m₀, P₀ = initial_moments(spec.init, p)
        a, c, s = _free_moments_affine(p, μ, m₀, P₀, y, free, k)
        rsetindex!(av, a, sidx + k - 1)
        rsetindex!(cv, c + s * rgetindex(x, index + k - 1), sidx + k - 1)
        ℓ += _maybe_logjac(flag, s)
    end
    flag isa TV.NoLogJac && return ℓ0, index + n, sidx + n
    return ℓ, index + n, sidx + n
end

function _color_chain_flat!(flag, y, x, index, sidx, av, cv, spec::IIDChainSpec, hp)
    ℓ = TV.logjac_zero(flag, eltype(x))
    dim = TV.dimension(spec.fnode)
    n = _nfree(spec)
    # Independent pointwise transforms; stage contiguously so the loop raises and the
    # scattered write is a single vectorized scatter (see `_color_chain_wrapped!`).
    ys = similar(y, n)
    @trace track_numbers = false for k in 1:n
        yi, ℓi, _ = TV.transform_with(flag, spec.fnode, x, index + (k - 1) * dim)
        rsetindex!(ys, yi, k)
        ℓ += ℓi
    end
    scatter_values!(y, spec.freeinds, ys)
    return ℓ, index + n * dim, sidx
end

@inline _whiten_specs_flat!(x, index, y, ::Tuple{}, hp) = index
@inline function _whiten_specs_flat!(x, index, y, specs::Tuple, hp)
    index = _whiten_chain_flat!(x, index, y, first(specs), hp)
    return _whiten_specs_flat!(x, index, y, Base.tail(specs), hp)
end

# Inverse of the centered wrapped lift. The coloring stores wrapped angles, so the
# coordinates are recovered by *unwrapping* the chain: each free point's coordinate is
# its own sheet anchor plus the wrapped step onto it, which places every point on the
# sheet the coloring's weight favours (steps within `(−π, π]`) instead of on whatever
# sheet the stored angle happened to land. Anchoring on the previous *coordinate* makes
# this a first-order recurrence, resolved with the same `affine_scan!` as the coloring
# rather than a loop-carried traced loop. `transform ∘ inverse` is the identity on
# wrapped angles; `inverse ∘ transform` is not (the sheet is not recoverable from the
# angle), as for the angle embedding.
# Inverse of the centered wrapped lift for a process whose transition mean is itself
# `2π`-periodic (`_periodic_mean`): `_cond_mean` then reads the previous point only through
# its *angle*, so the sheet the coloring's weight favours is pinned by the stored angles
# alone and the unwrap is pointwise — there is no recurrence to resolve, unlike the
# `WrappedBrownian` case below whose Φ ≡ 1 mean carries the previous point's sheet. The
# chain-opening point takes the centre of `_init_window`, the same window its weight used.
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

function _whiten_chain_flat!(x, index, y, spec::ChainUnit, hp)
    free = spec.free
    if is_wrapped(spec.process) && !spec.centered
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
    if spec.centered
        n = _nfree(spec)
        if is_wrapped(spec.process)
            # `_periodic_mean` is a trait on a concrete process type, so this is a
            # host-static branch.
            _periodic_mean(spec.process) &&
                return _whiten_chain_wrapped_pointwise!(x, index, y, spec, hp)
            return _whiten_chain_wrapped_centered!(x, index, y, spec)
        end
        @trace track_numbers = false for k in 1:n
            rsetindex!(x, rgetindex(y, rgetindex(free.tgt, k)), index + k - 1)
        end
        return index + n
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

function _color_chain_std!(y, x, index, sidx, av, cv, spec::ChainUnit, hp, space)
    # centered and wrapped chains are rejected at node construction (`_check_std_transportable`)
    vals = _unit_params(spec, hp)
    free = spec.free
    n = _nfree(spec)
    # Same affine split as the flat coloring (see `_color_chain_flat!`): coefficients
    # staged in a raisable loop; the caller's batched `affine_scan!` resolves them.
    @trace track_numbers = false for k in 1:n
        p = _unit_process(spec, vals, _hpidx(spec, free, k))
        μ = _process_mean(p)
        m₀, P₀ = initial_moments(spec.init, p)
        a, c, s = _free_moments_affine(p, μ, m₀, P₀, y, free, k)
        u = PT._clamp_unit(PT.space_cdf(space, rgetindex(x, index + k - 1)))
        rsetindex!(av, a, sidx + k - 1)
        rsetindex!(cv, c + s * PT.space_quantile(PT.StdNormal(), u), sidx + k - 1)
    end
    return index + n, sidx + n
end

function _color_chain_std!(y, x, index, sidx, av, cv, spec::IIDChainSpec, hp, space)
    tin = PT.transport_node(spec.dist, space)
    dim = PT.dimension(tin)
    n = _nfree(spec)
    # Independent pointwise transforms; stage contiguously so the loop raises and the
    # scattered write is a single vectorized scatter (see `_color_chain_wrapped!`).
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
        u = PT._clamp_unit(
            PT.space_cdf(PT.StdNormal(), (rgetindex(y, rgetindex(free.tgt, k)) - m) / s)
        )
        rsetindex!(x, PT.space_quantile(space, u), index + k - 1)
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

"""
    siteparams(x)

Extract the instrument parameter values from a single instrument-parameter sample. For
hierarchical samples `(params = SiteArray, hyperparams = NamedTuple)` (produced by
[`GaussMarkovSitePrior`](@ref) with fitted hyperparameters) this returns `x.params`;
for plain `SiteArray` samples it is the identity.
"""
@inline siteparams(x) = x
@inline siteparams(x::NamedTuple{(:params, :hyperparams)}) = x.params

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
# both to the data and between sites. (They used to be relative to each chain's own first
# point, which made `spec.ts` silently unmatchable against `Ti` and incomparable across
# sites.) Only *gaps* reach the process — `transition_moments`, the bridge, and the
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

function _chainspec(sp::GaussMarkovSitePrior, site, inds, ts, fixedpos)
    hpk = isempty(hyperprior(sp.process)) ? nothing : site
    return MarkovChainSpec(sp.process, sp.init, Val(hpk), inds, ts, fixedpos, sp.centered)
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

    T = _working_type(site_dists)
    finds, vals = reference_indices(array, smap, d.refant)
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
    # so its fitted fields are replicated for each site that uses it. Multi-frequency
    # chains of the same site share that site's hyperparameters.
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
