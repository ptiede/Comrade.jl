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

By default the prior is *whitened* (non-centered): each free point's latent coordinate is
an iid standard variate that is colored through the chain's exact conditional moments.
This removes the hyperparameter-gain funnel when hyperparameters are fitted and makes the
transport to the Std spaces exact, so both `asflat` and `ascube` work. Set
`centered = true` to instead use the gain values themselves as coordinates (the centered
parameterization), which can mix better when every point is strongly data-constrained;
centered priors support `asflat` only. Wrapped (circular) processes default to a third
parameterization: each free phase is embedded as two latent reals through the same angle
transform used by `DiagonalVonMises`, so the flat coordinates carry no `2π` ambiguity (a
phase wrap is a continuous winding of the latent point, never a jump) and HMC warmup
adaptation is robust even for weakly constrained phases. The embedding doubles the phase
dimension and its per-pair geometry is a rotated ellipse a diagonal mass matrix cannot
capture, so when every phase is strongly data-constrained `centered = true` mixes better:
it uses the raw angles as coordinates, at the cost of exactly `2π`-periodic flat
coordinates — a sampler sheet-hop on a weakly constrained point can destabilize warmup
variance adaptation. Wrapped chains support `asflat` only in either form.

This prior is intended for both gain amplitudes and phases. For phases the preferred
process is [`WrappedBrownian`](@ref) with `init = UniformInit()`: the prior is then
exactly `2π`-periodic, needs no separate circular offset term, and its diffusion
posterior is unbiased by phase wraps. A real-line process (e.g. `OrnsteinUhlenbeck` with
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

function GaussMarkovSitePrior(
        seg::TimeSegmentation, process::AbstractGaussMarkovProcess;
        init::AbstractInitialPrior = StationaryInit(), centered = false
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

@inline _selecthp(hp::NamedTuple, ::Val{nothing}) = (;)
@inline _selecthp(hp::NamedTuple, ::Val{K}) where {K} = getproperty(hp, K)

# ----- the chain log-density -------------------------------------------------

# Log-density of the chain's first point under the init marginal propagated a gap Δ0
# from the chain start (Δ0 = 0 for the full chain; the first fixed time for the
# conditioning subtraction in `chain_term`). For FixedInit the first point is always
# reference-fixed, so its (delta) term appears identically in ℓ(all) and ℓ(fixed subset)
# and is conventionally zero in both.
@inline function _first_point_term(init::AbstractInitialPrior, p, x1, Δ0)
    m1, P1 = marginal_moments(init, p, Δ0)
    return -(abs2(x1 - m1) / P1 + log(2π * P1)) / 2
end
@inline _first_point_term(::FixedInit, p, x1, Δ0) = zero(abs2(x1))
# Uniform on the circle, which is also invariant under the wrapped transitions, so the
# same constant serves the propagated subset marginal at any Δ0.
@inline _first_point_term(::UniformInit, p, x1, Δ0) = -log(oftype(x1, 2π))

function _gm_chain_logpdf(p::AbstractGaussMarkovProcess, init::AbstractInitialPrior, x, inds, ts, Δ0 = zero(eltype(ts)))
    x1 = rgetindex(x, rgetindex(inds, firstindex(inds)))
    ℓ = _first_point_term(init, p, x1, Δ0)
    μ = _process_mean(p)
    T2π = convert(eltype(x), 2π)
    @trace track_numbers = false for i in (firstindex(inds) + 1):lastindex(inds)
        xi = rgetindex(x, rgetindex(inds, i))
        xp = rgetindex(x, rgetindex(inds, i - 1))
        ℓ += _transition_logpdf(p, xi, xp, rgetindex(ts, i) - rgetindex(ts, i - 1), μ, T2π)
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
struct GaussMarkovChainDist{C <: NamedTuple, I <: AbstractVector{<:Integer}, F} <: Dists.ContinuousMultivariateDistribution
    chains::C
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

_scan_tgt(spec::MarkovChainSpec) =
    (is_wrapped(spec.process) || spec.centered) ? Int[] : spec.free.tgt
_scan_tgt(::IIDChainSpec) = Int[]

function GaussMarkovChainDist(chains::NamedTuple, fixedinds, fixedvals, len::Int)
    scantgt = reduce(vcat, map(_scan_tgt, collect(values(chains))); init = Int[])
    return GaussMarkovChainDist(chains, fixedinds, fixedvals, len, scantgt)
end

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
@inline (hp::_ChainFix)(spec::Union{MarkovChainSpec, IIDChainSpec}) = chain_term(spec, hp.x, hp.hp)

function _chain_logpdf(d::GaussMarkovChainDist, x::AbstractVector, hp::NamedTuple)
    fd = _ChainFix(x, hp)
    ls = map(fd, values(chainspecs(d)))
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
    return (; tgt, lidx, dtl, mskl, lfree, fidx, dtf, mskf, hasfix = nfix > 0)
end

_nfree(spec::MarkovChainSpec) = length(spec.free.tgt)
_nfree(spec::IIDChainSpec) = length(spec.freeinds)

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
    D = p.D
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

    step_walk(prev, Δt) = _wrap_angle(prev + sqrt(D * Δt) * randn(rng, T))

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
        w = _sample_winding(rng, Δθ, D * (tb - ta))
        e = x[inds[fa]] + Δθ + 2π * w  # unwrapped endpoint of this segment
        cur, tp = x[inds[fa]], ta
        for i in (fa + 1):(fb - 1)
            m = cur + (e - cur) * (ts[i] - tp) / (tb - tp)
            v = D * (ts[i] - tp) * (tb - ts[i]) / (tb - tp)
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
_flat_dim(spec::MarkovChainSpec) =
    (is_wrapped(spec.process) && !spec.centered) ? 2 * _nfree(spec) : _nfree(spec)
_flat_dim(spec::IIDChainSpec) = _nfree(spec) * TV.dimension(PT.transport_node(spec.dist, PT.TVFlat()))
_flat_dim(d::GaussMarkovChainDist) = sum(_flat_dim, values(chainspecs(d)); init = 0)

_std_dim(spec::MarkovChainSpec, space) = _nfree(spec)
_std_dim(spec::IIDChainSpec, space) = _nfree(spec) * PT.dimension(PT.transport_node(spec.dist, space))
_std_dim(d::GaussMarkovChainDist, space) = sum(Base.Fix2(_std_dim, space), values(chainspecs(d)); init = 0)

_std_transportable(spec::MarkovChainSpec) = !spec.centered && !is_wrapped(spec.process)
_std_transportable(::IIDChainSpec) = true
function _check_std_transportable(d::GaussMarkovChainDist)
    all(_std_transportable, values(chainspecs(d))) || throw(
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
    ℓ, index, _ = _color_specs_flat!(flag, y, x, index, 1, av, cv, values(chainspecs(d)), hp)
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
function _color_chain_wrapped!(flag, y, x, index, spec::MarkovChainSpec)
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

function _color_chain_flat!(flag, y, x, index, sidx, av, cv, spec::MarkovChainSpec, hp)
    free = spec.free
    n = _nfree(spec)
    ℓ0 = TV.logjac_zero(flag, eltype(x))
    n == 0 && return ℓ0, index, sidx
    if spec.centered
        # A straight copy of a contiguous flat segment to the scattered chain positions:
        # one vectorized scatter (a scattered-write loop would serialize under Reactant).
        scatter_values!(y, free.tgt, view(x, index:(index + n - 1)))
        return ℓ0, index + n, sidx
    end
    if is_wrapped(spec.process)
        ℓw, index = _color_chain_wrapped!(flag, y, x, index, spec)
        return ℓw, index, sidx
    end
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ = _process_mean(p)
    m₀, P₀ = initial_moments(spec.init, p)
    ℓ = zero(eltype(x))
    # The recurrence `y[tgt[k]] = m(y[tgt[k-1]]) + s*x[k]` is affine in the previous
    # free value, so split it: this loop computes the coefficients `(a, c)` — reading
    # only tables, hyperparameters, and fixed values, with contiguous writes, so it
    # raises under Reactant. The recurrence itself is resolved by the caller's single
    # batched `affine_scan!` over all chains.
    @trace track_numbers = false for k in 1:n
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

function _whiten_chain_flat!(x, index, y, spec::MarkovChainSpec, hp)
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
        @trace track_numbers = false for k in 1:n
            rsetindex!(x, rgetindex(y, rgetindex(free.tgt, k)), index + k - 1)
        end
        return index + n
    end
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ = _process_mean(p)
    m₀, P₀ = initial_moments(spec.init, p)
    @trace track_numbers = false for k in 1:_nfree(spec)
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
    index, _ = _color_specs_std!(y, x, index, 1, av, cv, values(chainspecs(d)), hp, space)
    _resolve_scan!(y, d, av, cv)
    return index
end

function _color_chain_std!(y, x, index, sidx, av, cv, spec::MarkovChainSpec, hp, space)
    # centered and wrapped chains are rejected at node construction (`_check_std_transportable`)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ = _process_mean(p)
    m₀, P₀ = initial_moments(spec.init, p)
    free = spec.free
    n = _nfree(spec)
    # Same affine split as the flat coloring (see `_color_chain_flat!`): coefficients
    # staged in a raisable loop; the caller's batched `affine_scan!` resolves them.
    @trace track_numbers = false for k in 1:n
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

function _whiten_chain_std!(x, index, y, spec::MarkovChainSpec, hp, space)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ = _process_mean(p)
    m₀, P₀ = initial_moments(spec.init, p)
    free = spec.free
    n = _nfree(spec)
    @trace track_numbers = false for k in 1:n
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
    return _whiten_specs_flat!(x, index, parent(y), values(chainspecs(t.dists)), (;))
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
    return _whiten_specs_std!(x, index, parent(y), values(chainspecs(t.dists)), (;), t.space)
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
    return _whiten_specs_flat!(x, index, parent(y.params), values(chainspecs(t.dists)), y.hyperparams)
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
    return _whiten_specs_std!(x, index, parent(y.params), values(chainspecs(t.dists)), y.hyperparams, t.space)
end

function PT.transport_node(d::ObservedHierarchicalArrayPrior, space::PT.AbstractStdDist)
    _check_std_transportable(d.dists)
    hnode = PT.transport_node(d.hyperprior, space)
    return StdWhitenedHierarchicalTransform(
        hnode, d.dists, d.sitemap, space, PT.dimension(hnode) + _std_dim(d.dists, space)
    )
end

# ----- construction from ArrayPrior --------------------------------------------

# Chain times in hours, relative to the first chain point (handles multi-day tracks via
# the mjd field).
function _chain_times(times, inds)
    t1 = times[first(inds)]
    ts = map(i -> (mjd(times[i]) - mjd(t1)) * 24.0 + (_center(times[i]) - _center(t1)), inds)
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
    chains = map(lookup(smap)) do inds
        s = smap.sites[first(inds)]
        sp = getproperty(site_dists, s)
        ts = _chain_times(smap.times, inds)
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
