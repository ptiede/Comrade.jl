export GaussMarkovChain

# ----- hyperparameter-aware transport shims -------------------------------------
#
# Chain transport nodes are built once (statically), but chains with fitted
# hyperparameters need the sampled hyperparameter values at *evaluation* time. These
# shims thread an `hp` NamedTuple through any transport node: the generic methods ignore
# it and delegate to the vanilla TransformVariables/ProbabilityTransports calls (which is
# correct for every wrapped `Distributions` node), while hyperparameter-aware nodes (the
# Gauss-Markov whitening nodes below) override them. The vanilla TV/PT entry points of
# those nodes delegate back to the shims with empty hyperparameters.

@inline _node_transform_with(flag::TV.LogJacFlag, node, x, index, hp) =
    TV.transform_with(flag, node, x, index)
@inline _node_inverse_at!(x, index, node, y, hp) = TV.inverse_at!(x, index, node, y)
@inline _node_pfwd(node, x, index, hp) = PT.pfwd_step(node, x, index)
@inline _node_pback!(x, index, node, y, hp) = PT.pback_step!(x, index, node, y)

# ----- the chain distribution ---------------------------------------------------

"""
    GaussMarkovChain(process::AbstractGaussMarkovProcess, times; centered = false, fixedpos = Int[], fixedvals = Float64[])

The distribution of a scalar Gauss-Markov `process` (e.g. [`OrnsteinUhlenbeck`](@ref),
[`Wiener`](@ref)) observed at the strictly increasing `times` (in hours), optionally
*exactly conditioned* on the values `fixedvals` at the positions `fixedpos` (indices into
`times`). This is an ordinary `Distributions.ContinuousMultivariateDistribution` over the
full chain vector (fixed positions included): the log-density is the density of the free
values given the fixed ones (`O(n)` via the process's exact discretization, using
`log p(free | fixed) = ℓ(all) − ℓ(fixed subset)`, valid since the restriction of a Markov
chain to a subset of times is Markov), `rand` samples the free values from the exact
conditional (bridge) law with the fixed values filled in, and the transport nodes
implement the exact whitening (non-centered) triangular map used by `asflat`/`ascube`.

Use [`condition`](@ref Comrade.condition) to add conditioning points to an existing
chain. Processes without a proper initial law (e.g. [`Wiener`](@ref)) require the first
point to be conditioned.

The `process` may contain `Distributions.Distribution` fields (fitted hyperparameters, see
[`AbstractGaussMarkovProcess`](@ref)); such a *templated* chain cannot be evaluated
directly — the instrument-prior machinery threads the sampled hyperparameter values
through it, or [`materialize`](@ref Comrade.materialize) substitutes them. When the
process is fully numeric the transition moments are precomputed at construction.

Set `centered = true` to make the latent coordinates the chain values themselves instead
of whitened standard variates; centered chains support the flat (`asflat`) transport only.
"""
struct GaussMarkovChain{P <: AbstractGaussMarkovProcess, T <: AbstractVector, F <: NamedTuple, S <: NamedTuple, C <: Union{Nothing, NamedTuple}} <: Dists.ContinuousMultivariateDistribution
    process::P            # process; fitted fields are Distribution placeholders
    ts::T                 # chain times in hours, shifted to start at zero
    centered::Bool        # latent coordinates are raw values instead of whitened variates
    dropinit::Bool        # drop the initial density term (first point fixed, or improper initial)
    fixedpos::Vector{Int} # conditioned positions, ascending
    fixedvals::Vector{Float64}
    free::F               # static per-free-point tables (see `_free_point_tables`)
    fsub::S               # (pos, ts) of the fixed subset
    tables::C             # precomputed transition moments when the process is numeric
end

function GaussMarkovChain(
        process::AbstractGaussMarkovProcess, times::AbstractVector{<:Real};
        centered::Bool = false,
        fixedpos::AbstractVector{<:Integer} = Int[],
        fixedvals::AbstractVector{<:Real} = Float64[]
    )
    n = length(times)
    n > 0 || throw(ArgumentError("GaussMarkovChain needs at least one time"))
    ts = float.(collect(times))
    ts .-= first(ts)
    for i in (firstindex(ts) + 1):lastindex(ts)
        ts[i] > ts[i - 1] || throw(ArgumentError("GaussMarkovChain times must be strictly increasing"))
    end
    length(fixedpos) == length(fixedvals) ||
        throw(ArgumentError("fixedpos and fixedvals must have equal length"))
    fp = collect(Int, fixedpos)
    fv = collect(float.(fixedvals))
    if !issorted(fp)
        p = sortperm(fp)
        fp, fv = fp[p], fv[p]
    end
    allunique(fp) || throw(ArgumentError("fixedpos must be unique"))
    isempty(fp) || (1 ≤ first(fp) && last(fp) ≤ n) ||
        throw(ArgumentError("fixedpos must be indices into times"))

    dropinit = !has_proper_initial(process) || (!isempty(fp) && first(fp) == 1)
    free = _free_point_tables(ts, fp)
    fsub = (pos = fp, ts = ts[fp])
    tables = isempty(hyperprior(process)) ? _transition_tables(process, ts, free, fsub, dropinit) : nothing
    return GaussMarkovChain(process, ts, centered, dropinit, fp, fv, free, fsub, tables)
end

Base.length(d::GaussMarkovChain) = length(d.ts)
Base.eltype(::GaussMarkovChain) = Float64
Dists.sampler(d::GaussMarkovChain) = d

_nfree(d::GaussMarkovChain) = length(d.free.tgt)

hyperprior(d::GaussMarkovChain) = hyperprior(d.process)

function materialize(d::GaussMarkovChain, hp::NamedTuple)
    isempty(hyperprior(d.process)) && return d
    return GaussMarkovChain(
        materialize(d.process, hp), d.ts;
        centered = d.centered, fixedpos = d.fixedpos, fixedvals = d.fixedvals
    )
end

function condition(d::GaussMarkovChain, pos::AbstractVector{<:Integer}, vals::AbstractVector{<:Real})
    # existing conditioning points win over newly added duplicates
    old = Set(d.fixedpos)
    keep = findall(!in(old), pos)
    return GaussMarkovChain(
        d.process, d.ts;
        centered = d.centered,
        fixedpos = vcat(d.fixedpos, collect(Int, pos[keep])),
        fixedvals = vcat(d.fixedvals, collect(float.(vals[keep])))
    )
end

@inline function _check_proper(d::GaussMarkovChain)
    has_proper_initial(d.process) && return nothing
    (!isempty(d.fixedpos) && first(d.fixedpos) == 1) && return nothing
    throw(
        ArgumentError(
            "$(nameof(typeof(d.process))) has no proper initial law, so the chain's first " *
                "point must be exactly conditioned. Use `anchored = true` in " *
                "GaussMarkovSitePrior, or condition the first point explicitly."
        )
    )
end

@inline function _check_numeric(d::GaussMarkovChain)
    isempty(hyperprior(d.process)) && return nothing
    throw(
        ArgumentError(
            "This GaussMarkovChain's process has fitted hyperparameters " *
                "($(keys(hyperprior(d.process)))). Substitute values with " *
                "`Comrade.materialize(d, hp)` before evaluating it directly."
        )
    )
end

# ----- static tables -------------------------------------------------------------

# Static per-free-point tables, computed once at construction: for free point `k` (in
# time order) the target position, the immediate left-neighbor position/gap/mask, and the
# next conditioned position/gap/mask. A missing neighbor is a `0.0` mask with a dummy
# self position and unit gap: the coefficient is *multiplied by the mask*, so the dummy
# read is scaled by an exact static zero. (An `Inf` gap would also zero the primal, but
# its τ-derivative is `0 * Inf = NaN` under AD.)
function _free_point_tables(ts, fixedpos)
    nfix = length(fixedpos)
    tgt = Int[]
    lidx = Int[]
    dtl = Float64[]
    mskl = Float64[]
    fidx = Int[]
    dtf = Float64[]
    mskf = Float64[]
    jfix = 1
    for i in eachindex(ts)
        while jfix ≤ nfix && fixedpos[jfix] < i
            jfix += 1
        end
        (jfix ≤ nfix && fixedpos[jfix] == i) && continue
        push!(tgt, i)
        if i > firstindex(ts)
            push!(lidx, i - 1)
            push!(dtl, ts[i] - ts[i - 1])
            push!(mskl, 1.0)
        else
            push!(lidx, i)
            push!(dtl, 1.0)
            push!(mskl, 0.0)
        end
        if jfix ≤ nfix
            b = fixedpos[jfix]
            push!(fidx, b)
            push!(dtf, ts[b] - ts[i])
            push!(mskf, 1.0)
        else
            push!(fidx, i)
            push!(dtf, 1.0)
            push!(mskf, 0.0)
        end
    end
    return (; tgt, lidx, dtl, mskl, fidx, dtf, mskf, hasfix = nfix > 0)
end

# Static transition-moment tables for fully numeric processes: a fixed-hyperparameter
# chain would otherwise recompute the same transcendentals (e.g. the OU Φ exponentials)
# on every log-density and transport evaluation. All tables are plain Float64 vectors
# read with the same gathers as the `free` index tables, so the traced loops are
# unchanged in structure. Masks are folded into the coefficients directly.
function _transition_tables(process::AbstractGaussMarkovProcess, ts, free, fsub, dropinit)
    p = materialize(process, (;))
    μ = process_mean(p)
    Φ1g, Q1g = transition_moments(p, one(eltype(ts)))
    (μ isa Number && Φ1g isa Number && Q1g isa Number) || throw(
        ArgumentError(
            "AbstractGaussMarkovProcess must have scalar state; vector-state processes " *
                "are not supported yet."
        )
    )
    if has_proper_initial(p)
        μ₀, P₀ = initial_moments(p)
        (μ₀ isa Number && P₀ isa Number) || throw(
            ArgumentError(
                "AbstractGaussMarkovProcess must have scalar state: initial_moments " *
                    "returned ($(typeof(μ₀)), $(typeof(P₀)))."
            )
        )
    else
        μ₀, P₀ = μ, 1.0 # dummies; never selected (dropinit, all free points have a left neighbor)
    end
    gaps(t) = [transition_moments(p, t[i] - t[i - 1]) for i in (firstindex(t) + 1):lastindex(t)]
    ΦQ = gaps(ts)
    fΦQ = gaps(fsub.ts)
    nf = length(free.tgt)
    Φ1 = Vector{Float64}(undef, nf)
    v1 = Vector{Float64}(undef, nf)
    a0 = Vector{Float64}(undef, nf)
    Φ2 = Vector{Float64}(undef, nf)
    Q2 = Vector{Float64}(undef, nf)
    for k in 1:nf
        Φl, Ql = transition_moments(p, free.dtl[k])
        Φ1[k] = free.mskl[k] * Φl
        v1[k] = free.mskl[k] * Ql + (1 - free.mskl[k]) * P₀
        a0[k] = (1 - free.mskl[k]) * (μ₀ - μ)
        Φf, Qf = transition_moments(p, free.dtf[k])
        Φ2[k] = free.mskf[k] * Φf
        Q2[k] = Qf
    end
    # marginal law at the first fixed time (for the ℓ(fixed) initial term), propagated
    # exactly from the initial law; only needed when the first point is free
    if !dropinit && !isempty(fsub.pos)
        Φd, Qd = transition_moments(p, first(fsub.ts) - first(ts))
        μf = μ + Φd * (μ₀ - μ)
        Pf = Φd^2 * P₀ + Qd
    else
        μf, Pf = μ, 1.0
    end
    return (;
        μ = float(μ), μ₀ = float(μ₀), P₀ = float(P₀), μf = float(μf), Pf = float(Pf),
        Φg = first.(ΦQ), Qg = last.(ΦQ),
        fΦg = first.(fΦQ), fQg = last.(fΦQ),
        Φ1, v1, a0, Φ2, Q2,
    )
end

EnzymeRules.inactive_type(::Type{<:GaussMarkovChain}) = true

# ----- moment sources ---------------------------------------------------------
#
# Where the chain log-density and walkers read the process moments and per-point
# conditional coefficients: chains with fitted hyperparameters materialize the process
# per evaluation (live), numeric chains read the precomputed `tables`. The choice
# dispatches on the tables field, so the (traced) loop bodies stay single-path.

struct LiveMoments{Pr, T}
    p::Pr
    μ::T
    μ₀::T
    P₀::T
end
struct TableMoments{T <: NamedTuple}
    t::T
end

@inline _moments(d::GaussMarkovChain, hp::NamedTuple) = _moments_impl(d.tables, d, hp)
@inline _moments_impl(t::NamedTuple, d, hp) = TableMoments(t)
@inline function _moments_impl(::Nothing, d, hp)
    p = materialize(d.process, hp)
    μ = process_mean(p)
    if has_proper_initial(p)
        μ₀, P₀ = initial_moments(p)
    else
        μ₀, P₀ = μ, one(μ)
    end
    μp, μ₀p, P₀p = promote(μ, μ₀, P₀)
    return LiveMoments(p, μp, μ₀p, P₀p)
end

@inline _pmean(m::LiveMoments) = m.μ
@inline _pmean(m::TableMoments) = m.t.μ
@inline _init_moments(m::LiveMoments) = (m.μ₀, m.P₀)
@inline _init_moments(m::TableMoments) = (m.t.μ₀, m.t.P₀)

# marginal law at the first fixed time, propagated from the initial law
@inline function _finit_moments(m::LiveMoments, d::GaussMarkovChain)
    Φd, Qd = transition_moments(m.p, first(d.fsub.ts) - first(d.ts))
    return m.μ + Φd * (m.μ₀ - m.μ), Φd^2 * m.P₀ + Qd
end
@inline _finit_moments(m::TableMoments, d::GaussMarkovChain) = (m.t.μf, m.t.Pf)

# masked conditional coefficients of free point `k` (see `_free_point_tables`)
@inline function _left_coeffs(m::LiveMoments, free, k)
    mskl = rgetindex(free.mskl, k)
    Φl, Ql = transition_moments(m.p, rgetindex(free.dtl, k))
    return mskl * Φl, mskl * Ql + (1 - mskl) * m.P₀, (1 - mskl) * (m.μ₀ - m.μ)
end
@inline _left_coeffs(m::TableMoments, free, k) =
    (rgetindex(m.t.Φ1, k), rgetindex(m.t.v1, k), rgetindex(m.t.a0, k))
@inline function _right_coeffs(m::LiveMoments, free, k)
    Φf, Qf = transition_moments(m.p, rgetindex(free.dtf, k))
    return rgetindex(free.mskf, k) * Φf, Qf
end
@inline _right_coeffs(m::TableMoments, free, k) =
    (rgetindex(m.t.Φ2, k), rgetindex(m.t.Q2, k))

# Per-gap transition moments (gap `i` sits between points `i-1` and `i`), either computed
# live or gathered from the static tables. `_gaps`/`_fsubgaps` select the main-chain vs
# fixed-subset table pair.
struct LiveGaps{Pr, T}
    p::Pr
    ts::T
end
struct TableGaps{V}
    Φ::V
    Q::V
end

@inline _gapΦQ(g::LiveGaps, i) =
    transition_moments(g.p, rgetindex(g.ts, i) - rgetindex(g.ts, i - 1))
@inline _gapΦQ(g::TableGaps, i) = (rgetindex(g.Φ, i - 1), rgetindex(g.Q, i - 1))

@inline _gaps(m::LiveMoments, d::GaussMarkovChain) = LiveGaps(m.p, d.ts)
@inline _gaps(m::TableMoments, d::GaussMarkovChain) = TableGaps(m.t.Φg, m.t.Qg)
@inline _fsubgaps(m::LiveMoments, d::GaussMarkovChain) = LiveGaps(m.p, d.fsub.ts)
@inline _fsubgaps(m::TableMoments, d::GaussMarkovChain) = TableGaps(m.t.fΦg, m.t.fQg)

# ----- the chain log-density -------------------------------------------------

@inline _gauss_loglik(μ, P, x) = -(abs2(x - μ) / P + log(2π * P)) / 2

# Transition terms of the whole chain: O(n) and Enzyme-safe (sequential loop). The
# previous point is *re-read* each iteration rather than carried across iterations: the
# only loop-carried state is then the accumulator, which lets Reactant raise the `@trace`
# loop to fused vector ops instead of a serialized `stablehlo.while` (a carried value
# defeats the raising pass). On the CPU the extra read is free. `track_numbers = false`
# keeps `@trace` from promoting captured literals (e.g. `π`). The gap source is
# dispatch-static per chain, so the loop body stays single-path.
function _gm_transitions(μ, gaps, x)
    ℓ = zero(μ * rgetindex(x, firstindex(x)))
    T2π = convert(eltype(x), 2π)
    @trace track_numbers = false for i in (firstindex(x) + 1):lastindex(x)
        Φ, Q = _gapΦQ(gaps, i)
        xi = rgetindex(x, i)
        xp = rgetindex(x, i - 1)
        ℓ -= (abs2(xi - μ - Φ * (xp - μ)) / Q + log(T2π * Q)) / 2
    end
    return ℓ
end

# same, over the subset of points `pos` (the fixed subset)
function _gm_transitions(μ, gaps, x, pos)
    ℓ = zero(μ * rgetindex(x, rgetindex(pos, firstindex(pos))))
    T2π = convert(eltype(x), 2π)
    @trace track_numbers = false for i in (firstindex(pos) + 1):lastindex(pos)
        Φ, Q = _gapΦQ(gaps, i)
        xi = rgetindex(x, rgetindex(pos, i))
        xp = rgetindex(x, rgetindex(pos, i - 1))
        ℓ -= (abs2(xi - μ - Φ * (xp - μ)) / Q + log(T2π * Q)) / 2
    end
    return ℓ
end

# Exact conditioning on the fixed values: the restriction of an order-1 Markov process to
# a subset of times is Markov with the same transition law over the larger gaps, so
# log p(free | fixed) = ℓ(all) − ℓ(fixed subset). The subtracted term depends on the
# hyperparameters even though the fixed values are constants — dropping it would bias the
# hyperparameter posterior. When the first point is fixed (or the initial law is
# improper) the initial density terms are dropped from both: they cancel exactly.
function _cd_logpdf(d::GaussMarkovChain, x::AbstractVector, hp::NamedTuple)
    _check_proper(d)
    m = _moments(d, hp)
    μ = _pmean(m)
    ℓ = _gm_transitions(μ, _gaps(m, d), x)
    if !d.dropinit
        μ₀, P₀ = _init_moments(m)
        ℓ += _gauss_loglik(μ₀, P₀, rgetindex(x, firstindex(x)))
    end
    isempty(d.fixedpos) && return ℓ
    ℓf = _gm_transitions(μ, _fsubgaps(m, d), x, d.fsub.pos)
    if !d.dropinit
        μf, Pf = _finit_moments(m, d)
        ℓf += _gauss_loglik(μf, Pf, rgetindex(x, first(d.fsub.pos)))
    end
    return ℓ - ℓf
end

function Dists.logpdf(d::GaussMarkovChain, x::AbstractVector{<:Number})
    _check_numeric(d)
    return _cd_logpdf(d, x, (;))
end

# ----- the conditional (bridge) moments of a free point --------------------------
#
# Every per-point operation on a chain — the bridge sampler, the whitened coloring, its
# inverse, and the Std-space transports — walks the *free points only* and needs the same
# exact conditional moments: at free point `k`, condition on the previously realized
# value (or the initial law when there is none) and the *next conditioned* value, which
# by the Markov property is the exact conditional. In centered variables the left
# prediction is `N(a, v₁)` and the right neighbor contributes a Gaussian observation
# `y₂ ~ N(Φ₂ y, Q₂)`, so the posterior is the standard precision-weighted combination
#     1/v = 1/v₁ + Φ₂²/Q₂,   m = v (a/v₁ + Φ₂ y₂ / Q₂),
# with masked coefficients when a neighbor is missing (recovering the one-sided and
# initial-law conditionals). No stationarity is assumed. The walk is branchless over the
# static `free` tables so the same loop serves the CPU and, via `@trace`, compiles to a
# single while loop under Reactant. The `hasfix` branch is static (per chain), so chains
# without conditioned points skip the Φ₂ transition entirely — halving the transcendental
# cost for refant-free priors.
@inline function _free_moments(m, y, free, k)
    μ = _pmean(m)
    Φ₁, v₁, a0 = _left_coeffs(m, free, k)
    a = a0 + Φ₁ * (rgetindex(y, rgetindex(free.lidx, k)) - μ)
    if free.hasfix
        Φ₂, Q₂ = _right_coeffs(m, free, k)
        y₂ = rgetindex(y, rgetindex(free.fidx, k)) - μ
        w = Φ₂ / Q₂
        v = inv(inv(v₁) + Φ₂ * w)
        return μ + v * (a / v₁ + w * y₂), sqrt(v)
    end
    return μ + a, sqrt(v₁)
end

# Zero-initialize the chain vector and scatter the conditioned values. The zero fill
# matters: the branchless walk reads masked dummy neighbor positions, which must not be
# uninitialized memory (0 * NaN = NaN).
function _fill_fixed_chain!(y, d::GaussMarkovChain)
    fill!(y, zero(eltype(y)))
    @inbounds for (i, v) in zip(d.fixedpos, d.fixedvals)
        rsetindex!(y, v, i)
    end
    return y
end

# ----- exact conditional sampling (bridge walk) -----------------------------------
# rand-only code, never in the logpdf hot path. `_fill_fixed_chain!` zero-initializes,
# which the masked dummy reads of `_free_moments` require.

function _cd_rand!(rng::Random.AbstractRNG, d::GaussMarkovChain, x::AbstractVector, hp::NamedTuple)
    _check_proper(d)
    _fill_fixed_chain!(x, d)
    m = _moments(d, hp)
    free = d.free
    @inbounds for k in 1:_nfree(d)
        mm, s = _free_moments(m, x, free, k)
        x[free.tgt[k]] = mm + s * randn(rng, typeof(mm))
    end
    return x
end

function Dists._rand!(rng::Random.AbstractRNG, d::GaussMarkovChain, x::AbstractVector{<:Real})
    _check_numeric(d)
    return _cd_rand!(rng, d, x, (;))
end

# ----- flat and Std transport nodes: whitening / coloring ----------------------
#
# A whitened (non-centered) chain's latent coordinates are iid standard variates that are
# *colored* through the triangular map above: each free point is `m + s * z`. The map is
# exact, so it transports the Std spaces (`ascube`/`StdNormal`) with no density
# bookkeeping, while in flat space TV threads the `Σ log s` Jacobian and the
# constrained-space chain logpdf is unchanged. Centered chains keep their raw values as
# coordinates (identity, zero Jacobian) and therefore only support the flat space, where
# the target logpdf is evaluated.

# --- latent-space adapters: the color/whiten walkers differ between the flat
# (TransformVariables) and Std (ProbabilityTransports) transports only in (a) how a
# latent coordinate is converted to/from a standard variate and (b) whether a
# log-Jacobian is threaded. The adapters are zero-size and dispatch-static, so the shared
# (traced) loop bodies stay single-path.

struct FlatLatent{L <: TV.LogJacFlag}
    flag::L
end
struct StdLatent{S <: PT.AbstractStdDist}
    space::S
end

# read latent coordinate `x[index]` as a standard normal variate
@inline _pull_std(::FlatLatent, x, index) = rgetindex(x, index)
@inline function _pull_std(a::StdLatent, x, index)
    u = PT._clamp_unit(PT.space_cdf(a.space, rgetindex(x, index)))
    return PT.space_quantile(PT.StdNormal(), u)
end
# write standard normal variate `z` back to a latent coordinate
@inline _push_std(::FlatLatent, z) = z
@inline function _push_std(a::StdLatent, z)
    u = PT._clamp_unit(PT.space_cdf(PT.StdNormal(), z))
    return PT.space_quantile(a.space, u)
end

# Skip the `log` on Jacobian-free paths; dispatch is on a compile-time flag, so the
# traced loop body stays single-path. Std transports are exact (no Jacobian).
@inline _maybe_logjac(::TV.NoLogJac, s) = zero(s)
@inline _maybe_logjac(::TV.LogJac, s) = log(s)
@inline _maybe_logjac(a::FlatLatent, s) = _maybe_logjac(a.flag, s)
@inline _maybe_logjac(::StdLatent, s) = zero(s)

@inline _logjac_zero(a::FlatLatent, ::Type{T}) where {T} = TV.logjac_zero(a.flag, T)
@inline _logjac_zero(::StdLatent, ::Type{T}) where {T} = zero(T)

@inline _use_logjac(a::FlatLatent) = !(a.flag isa TV.NoLogJac)
@inline _use_logjac(::StdLatent) = false

# --- the chain walkers ---

function _color_chain!(a, y, x, index, d::GaussMarkovChain, hp)
    free = d.free
    n = _nfree(d)
    ℓ0 = _logjac_zero(a, eltype(x))
    n == 0 && return ℓ0, index
    if d.centered
        # centered chains are flat-only (`transport_node` refuses the Std spaces)
        @trace track_numbers = false for k in 1:n
            rsetindex!(y, rgetindex(x, index + k - 1), rgetindex(free.tgt, k))
        end
        return ℓ0, index + n
    end
    m = _moments(d, hp)
    ℓ = zero(eltype(x))
    @trace track_numbers = false for k in 1:n
        mm, s = _free_moments(m, y, free, k)
        rsetindex!(y, mm + s * _pull_std(a, x, index + k - 1), rgetindex(free.tgt, k))
        ℓ += _maybe_logjac(a, s)
    end
    _use_logjac(a) || return ℓ0, index + n
    return ℓ, index + n
end

function _whiten_chain!(a, x, index, y, d::GaussMarkovChain, hp)
    free = d.free
    if d.centered
        @inbounds for k in 1:_nfree(d)
            x[index + k - 1] = y[free.tgt[k]]
        end
        return index + _nfree(d)
    end
    m = _moments(d, hp)
    @inbounds for k in 1:_nfree(d)
        mm, s = _free_moments(m, y, free, k)
        x[index + k - 1] = _push_std(a, (y[free.tgt[k]] - mm) / s)
    end
    return index + _nfree(d)
end

# --- the transport nodes ---

struct GaussMarkovWhitenTransform{D <: GaussMarkovChain} <: TV.VectorTransform
    d::D
end

TV.dimension(t::GaussMarkovWhitenTransform) = _nfree(t.d)
TV.inverse_eltype(::GaussMarkovWhitenTransform, x::Type) = eltype(x)

TV.transform_with(flag::TV.LogJacFlag, t::GaussMarkovWhitenTransform, x, index) =
    _node_transform_with(flag, t, x, index, (;))
TV.inverse_at!(x::AbstractArray, index, t::GaussMarkovWhitenTransform, y::AbstractVector) =
    _node_inverse_at!(x, index, t, y, (;))

function _node_transform_with(flag::TV.LogJacFlag, t::GaussMarkovWhitenTransform, x, index, hp)
    d = t.d
    y = similar(x, length(d))
    _fill_fixed_chain!(y, d)
    return _color_chain_ret(FlatLatent(flag), y, x, index, d, hp)
end

@inline function _color_chain_ret(a, y, x, index, d, hp)
    ℓ, index = _color_chain!(a, y, x, index, d, hp)
    return y, ℓ, index
end

function _node_inverse_at!(x, index, t::GaussMarkovWhitenTransform, y, hp)
    return _whiten_chain!(FlatLatent(TV.NoLogJac()), x, index, y, t.d, hp)
end

struct StdGaussMarkovWhitenTransform{D <: GaussMarkovChain, S <: PT.AbstractStdDist} <: PT.AbstractTransport
    d::D
    space::S
end

PT.dimension(t::StdGaussMarkovWhitenTransform) = _nfree(t.d)
PT.pback_eltype(::StdGaussMarkovWhitenTransform) = Float64

PT.pfwd_step(t::StdGaussMarkovWhitenTransform, x, index) = _node_pfwd(t, x, index, (;))
PT.pback_step!(x::AbstractVector, index, t::StdGaussMarkovWhitenTransform, y) =
    _node_pback!(x, index, t, y, (;))

function _node_pfwd(t::StdGaussMarkovWhitenTransform, x, index, hp)
    d = t.d
    y = similar(x, float(eltype(x)), length(d))
    _fill_fixed_chain!(y, d)
    _, index = _color_chain!(StdLatent(t.space), y, x, index, d, hp)
    return y, index
end

function _node_pback!(x, index, t::StdGaussMarkovWhitenTransform, y, hp)
    return _whiten_chain!(StdLatent(t.space), x, index, y, t.d, hp)
end

EnzymeRules.inactive_type(::Type{<:GaussMarkovWhitenTransform}) = true
EnzymeRules.inactive_type(::Type{<:StdGaussMarkovWhitenTransform}) = true

function PT.transport_node(d::GaussMarkovChain, ::PT.TVFlat)
    _check_proper(d)
    return GaussMarkovWhitenTransform(d)
end

function PT.transport_node(d::GaussMarkovChain, space::PT.AbstractStdDist)
    _check_proper(d)
    d.centered && throw(
        ArgumentError(
            "GaussMarkovChain(...; centered = true) cannot be transported to the Std " *
                "spaces (ascube/StdNormal): the centered parameterization needs the target " *
                "log-density, which those transports never evaluate. Use the default " *
                "whitened form (centered = false) or asflat()."
        )
    )
    return StdGaussMarkovWhitenTransform(d, space)
end
