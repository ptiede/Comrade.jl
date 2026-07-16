export GaussMarkovSitePrior

"""
    GaussMarkovSitePrior(seg::TimeSegmentation, process::AbstractGaussMarkovProcess; centered = false, anchored = false)

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
centered priors support `asflat` only.

This prior is intended for both gain amplitudes and phases. Phases are modeled as
*unwrapped* reals (the likelihood only uses `exp(iθ)`, so wrapping is a display concern),
and the `phase=true` cumulative reparameterization of `ArrayPrior` is superseded by this
prior: `phase` must be `false`. Reference stations (`refant`) are handled by exact
conditioning of the chain on the fixed values, which works with scattered fixed indices
such as those produced by `SEFDReference`.

Setting `anchored = true` additionally conditions each site's chain to be exactly zero at
its first time, so the prior describes the *evolution* of a quantity relative to its
starting value (the continuous-time analogue of the cumulative `phase=true`
parameterization). This is important when the chain models an unwrapped phase fluctuation
on top of a separate offset term: without anchoring the chain's overall level trades off
against the offset, and the likelihood's `2π` periodicity turns that redundancy into
spurious posterior modes at `2π`-shifted levels. Anchoring removes the level freedom
entirely, so the wrap ambiguity lives only in the (circular) offset while genuine
continuous drift is still expressed by the chain.

# Example

```julia
lg ~ ArrayPrior(GaussMarkovSitePrior(IntegSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0))))
gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 3.0, τ = 1.0)); refant = SEFDReference(0.0))
## phase-fluctuation term meant to be combined with a circular per-track offset
dgp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = 1.0); anchored = true); refant = SEFDReference(0.0))
```
"""
struct GaussMarkovSitePrior{S <: TimeSegmentation, P <: AbstractGaussMarkovProcess} <: AbstractSitePrior
    seg::S
    process::P
    centered::Bool
    anchored::Bool
end

function GaussMarkovSitePrior(seg::TimeSegmentation, process::AbstractGaussMarkovProcess; centered = false, anchored = false)
    return GaussMarkovSitePrior(seg, process, centered, anchored)
end

# ----- per-(site, freq) chain specifications --------------------------------
#
# Built once at `set_array` time; everything in a spec is constant during sampling
# (indices, times, and the process *template* whose fitted fields are placeholders).

struct MarkovChainSpec{P <: AbstractGaussMarkovProcess, K, V <: AbstractVector{<:Integer}, T <: AbstractVector, F <: NamedTuple, G <: NamedTuple}
    process::P     # template; fitted fields are Distribution placeholders
    hpsel::Val{K}  # the site whose hyperparameters this chain reads, or nothing if fully fixed
    inds::V        # flat indices into the full parameter vector, ascending in time
    ts::T          # chain times in hours (strictly increasing)
    centered::Bool # coordinates are raw values (centered) or whitened standard variates
    free::F        # static per-free-point tables (see `_free_point_tables`)
    fsub::G        # inds/ts restricted to the fixed points, materialized (views do not trace)
end

function MarkovChainSpec(
        process::AbstractGaussMarkovProcess, hpsel::Val, inds, ts, fixedpos, centered::Bool
    )
    return MarkovChainSpec(
        process, hpsel, inds, ts, centered,
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

# Stationary start plus exact transition terms: O(n) and Enzyme-safe (sequential loop,
# `site_sum` precedent). The previous point is *re-read* each iteration rather than
# carried across iterations: the only loop-carried state is then the accumulator, which
# lets Reactant raise the `@trace` loop to fused vector ops instead of a serialized
# `stablehlo.while` (a carried value defeats the raising pass; see the Reactant issue
# MRE). On the CPU the extra read is free. `track_numbers = false` keeps `@trace` from
# promoting captured literals (e.g. `π`).
function _gm_chain_logpdf(p::AbstractGaussMarkovProcess, x, inds, ts)
    μ, P = stationary_moments(p)
    ℓ = -(abs2(rgetindex(x, rgetindex(inds, firstindex(inds))) - μ) / P + log(2π * P)) / 2
    T2π = convert(eltype(x), 2π)
    @trace track_numbers = false for i in (firstindex(inds) + 1):lastindex(inds)
        Φ, Q = transition_moments(p, rgetindex(ts, i) - rgetindex(ts, i - 1))
        xi = rgetindex(x, rgetindex(inds, i))
        xp = rgetindex(x, rgetindex(inds, i - 1))
        ℓ -= (abs2(xi - μ - Φ * (xp - μ)) / Q + log(T2π * Q)) / 2
    end
    return ℓ
end

# Exact conditioning on the reference-fixed values: the restriction of an order-1 Markov
# process to a subset of times is Markov with the same transition law over the larger
# gaps, so log p(free | fixed) = ℓ(all) − ℓ(fixed subset). The subtracted term depends on
# the hyperparameters even though the fixed values are constants — dropping it would bias
# the hyperparameter posterior.
@inline function chain_term(spec::MarkovChainSpec, x, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    ℓ = _gm_chain_logpdf(p, x, spec.inds, spec.ts)
    isempty(spec.fsub.inds) && return ℓ
    return ℓ - _gm_chain_logpdf(p, x, spec.fsub.inds, spec.fsub.ts)
end

@inline function chain_term(spec::IIDChainSpec, x, hp)
    B = Base.Fix1(Dists.logpdf, spec.dist)
    G = Base.Fix1(rgetindex, x)
    F = B ∘ G
    # `init`: `freeinds` is empty when every point of the site is reference-fixed
    return sum(F, spec.freeinds; init = zero(eltype(x)))
end

# ----- the full-vector distribution ------------------------------------------

"""
    GaussMarkovChainDist

The observed distribution over the full instrument parameter vector implied by
[`GaussMarkovSitePrior`](@ref): a set of independent per-(site, frequency) Gauss-Markov
chains (plus IID terms for `IIDSitePrior` overrides), exactly conditioned on the
reference-fixed indices. This is an internal type constructed by `ObservedArrayPrior`.
"""
struct GaussMarkovChainDist{C <: NamedTuple, I <: AbstractVector{<:Integer}, F} <: Dists.ContinuousMultivariateDistribution
    chains::C
    fixedinds::I
    fixedvals::F
    len::Int
end

Base.length(d::GaussMarkovChainDist) = d.len
Base.eltype(::GaussMarkovChainDist) = Float64
Dists.sampler(d::GaussMarkovChainDist) = d

chainspecs(d::GaussMarkovChainDist) = d.chains
EnzymeRules.inactive(::typeof(chainspecs), args...) = nothing

function _chain_logpdf(d::GaussMarkovChainDist, x::AbstractVector, hp::NamedTuple)
    ls = map(spec -> chain_term(spec, x, hp), chainspecs(d))
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
# exact conditional moments: at free point `k`, condition on the previously realized
# value (gap Δ₁ ⇒ Φ₁) and the *next reference-fixed* value (gap Δ₂ ⇒ Φ₂), which by the
# Markov property is the exact conditional. In centered variables y = x − μ:
#     mean = [Φ₁(1−Φ₂²) y_prev + Φ₂(1−Φ₁²) y_next] / (1 − Φ₁²Φ₂²)
#     var  = P (1−Φ₁²)(1−Φ₂²) / (1 − Φ₁²Φ₂²)
# with Φ = 0 when the neighbor is missing, which recovers the stationary and one-sided
# conditionals. The walk is branchless over the static `spec.free` tables so the same
# loop serves the CPU and, via `@trace`, compiles to a single while loop under Reactant
# (the fixed-point skip logic of a sequential walk cannot be loop-carried in a traced
# while).

@inline function _bridge_moments(μ, P, Φ₁, y₁, Φ₂, y₂)
    den = 1 - Φ₁^2 * Φ₂^2
    m = μ + (Φ₁ * (1 - Φ₂^2) * y₁ + Φ₂ * (1 - Φ₁^2) * y₂) / den
    s = sqrt(P * (1 - Φ₁^2) * (1 - Φ₂^2) / den)
    return m, s
end

# Conditional moments of free point `k`, reading the realized values from `y`. The
# `hasfix` branch is static (per chain), so chains without reference-fixed points skip
# the Φ₂ transition entirely — halving the transcendental cost for refant-free priors —
# while still compiling to a single traced loop.
@inline function _free_moments(p, μ, P, y, free, k)
    Φ₁ = rgetindex(free.mskl, k) * first(transition_moments(p, rgetindex(free.dtl, k)))
    y₁ = rgetindex(y, rgetindex(free.lidx, k)) - μ
    if free.hasfix
        Φ₂ = rgetindex(free.mskf, k) * first(transition_moments(p, rgetindex(free.dtf, k)))
        y₂ = rgetindex(y, rgetindex(free.fidx, k)) - μ
        return _bridge_moments(μ, P, Φ₁, y₁, Φ₂, y₂)
    end
    return _bridge_moments(μ, P, Φ₁, y₁, zero(Φ₁), zero(y₁))
end

# Zero-initialize the full parameter vector and scatter the reference-fixed values. The
# zero fill matters: the branchless walk reads masked dummy neighbor positions, which
# must not be uninitialized memory (0 * NaN = NaN).
function _fill_fixed!(y, d::GaussMarkovChainDist)
    fill!(y, zero(eltype(y)))
    @inbounds for (i, v) in zip(d.fixedinds, d.fixedvals)
        rsetindex!(y, v, i)
    end
    return y
end

# Static per-free-point tables for a Markov chain, computed once at construction: for
# free point `k` (in coordinate order) the target flat index, the left-neighbor flat
# index/gap/mask, and the next-reference-fixed flat index/gap/mask. A missing neighbor
# is a `0.0` mask with a dummy self index and unit gap: the transition coefficient is
# *multiplied by the mask*, so the dummy read is scaled by an exact static zero. (An
# `Inf` gap would also zero the primal, but its τ-derivative is `0 * Inf = NaN` under
# AD.)
function _free_point_tables(inds, ts, fixedpos)
    nfix = length(fixedpos)
    tgt = Int[]
    lidx = Int[]
    dtl = Float64[]
    mskl = Float64[]
    fidx = Int[]
    dtf = Float64[]
    mskf = Float64[]
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
            push!(mskl, 1.0)
        else
            push!(lidx, inds[i])
            push!(dtl, 1.0)
            push!(mskl, 0.0)
        end
        if jfix ≤ nfix
            b = fixedpos[jfix]
            push!(fidx, inds[b])
            push!(dtf, ts[b] - ts[i])
            push!(mskf, 1.0)
        else
            push!(fidx, inds[i])
            push!(dtf, 1.0)
            push!(mskf, 0.0)
        end
    end
    return (; tgt, lidx, dtl, mskl, fidx, dtf, mskf, hasfix = nfix > 0)
end

_nfree(spec::MarkovChainSpec) = length(spec.free.tgt)
_nfree(spec::IIDChainSpec) = length(spec.freeinds)

# ----- exact conditional sampling (OU bridge) -----------------------------------
# rand-only code, never in the logpdf hot path. `_fill_fixed!` zero-initializes, which
# the masked dummy reads of `_free_moments` require.

function _rand_chains!(rng::Random.AbstractRNG, x::AbstractVector, d::GaussMarkovChainDist, hp::NamedTuple)
    _fill_fixed!(x, d)
    foreach(spec -> _rand_chain!(rng, x, spec, hp), chainspecs(d))
    return x
end

function _rand_chain!(rng::Random.AbstractRNG, x, spec::IIDChainSpec, hp)
    @inbounds for i in spec.freeinds
        x[i] = rand(rng, spec.dist)
    end
    return x
end

function _rand_chain!(rng::Random.AbstractRNG, x, spec::MarkovChainSpec, hp)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ, P = stationary_moments(p)
    free = spec.free
    @inbounds for k in 1:_nfree(spec)
        m, s = _free_moments(p, μ, P, x, free, k)
        x[free.tgt[k]] = m + s * randn(rng, typeof(m))
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

_flat_dim(spec::MarkovChainSpec) = _nfree(spec)
_flat_dim(spec::IIDChainSpec) = _nfree(spec) * TV.dimension(PT.transport_node(spec.dist, PT.TVFlat()))
_flat_dim(d::GaussMarkovChainDist) = sum(_flat_dim, values(chainspecs(d)); init = 0)

_std_dim(spec::MarkovChainSpec, space) = _nfree(spec)
_std_dim(spec::IIDChainSpec, space) = _nfree(spec) * PT.dimension(PT.transport_node(spec.dist, space))
_std_dim(d::GaussMarkovChainDist, space) = sum(Base.Fix2(_std_dim, space), values(chainspecs(d)); init = 0)

_std_transportable(spec::MarkovChainSpec) = !spec.centered
_std_transportable(::IIDChainSpec) = true
function _check_std_transportable(d::GaussMarkovChainDist)
    all(_std_transportable, values(chainspecs(d))) || throw(
        ArgumentError(
            "GaussMarkovSitePrior(...; centered = true) cannot be transported to the Std " *
                "spaces (ascube/StdNormal): the centered parameterization needs the target " *
                "log-density, which those transports never evaluate. Use the default " *
                "whitened form (centered = false) or asflat()."
        )
    )
    return nothing
end

# --- flat (TransformVariables) chain walkers: thread (logjac, index) over the
# heterogeneous spec tuple with `Base.tail` recursion (TV tuple-transform precedent).

@inline _color_specs_flat!(flag, y, x, index, ::Tuple{}, hp) = TV.logjac_zero(flag, eltype(x)), index
@inline function _color_specs_flat!(flag, y, x, index, specs::Tuple, hp)
    ℓ1, index = _color_chain_flat!(flag, y, x, index, first(specs), hp)
    ℓ2, index = _color_specs_flat!(flag, y, x, index, Base.tail(specs), hp)
    return ℓ1 + ℓ2, index
end

# Skip the `log` on the value-only path; dispatch is on the compile-time flag, so the
# traced loop body stays single-path.
@inline _maybe_logjac(::TV.NoLogJac, s) = zero(s)
@inline _maybe_logjac(::TV.LogJac, s) = log(s)

function _color_chain_flat!(flag, y, x, index, spec::MarkovChainSpec, hp)
    free = spec.free
    n = _nfree(spec)
    ℓ0 = TV.logjac_zero(flag, eltype(x))
    n == 0 && return ℓ0, index
    if spec.centered
        y[free.tgt] = @view(x[index:(index + n - 1)])
        # rsetindex!(y, rgetindex(x, index:(index + n - 1)), free.tgt)
        return ℓ0, index + n
    end
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ, P = stationary_moments(p)
    ℓ = zero(eltype(x))
    @trace track_numbers = false for k in 1:n
        m, s = _free_moments(p, μ, P, y, free, k)
        rsetindex!(y, m + s * rgetindex(x, index + k - 1), rgetindex(free.tgt, k))
        ℓ += _maybe_logjac(flag, s)
    end
    flag isa TV.NoLogJac && return ℓ0, index + n
    return ℓ, index + n
end

function _color_chain_flat!(flag, y, x, index, spec::IIDChainSpec, hp)
    ℓ = TV.logjac_zero(flag, eltype(x))
    @inbounds for i in spec.freeinds
        yi, ℓi, index = TV.transform_with(flag, spec.fnode, x, index)
        rsetindex!(y, yi, i)
        ℓ += ℓi
    end
    return ℓ, index
end

@inline _whiten_specs_flat!(x, index, y, ::Tuple{}, hp) = index
@inline function _whiten_specs_flat!(x, index, y, specs::Tuple, hp)
    index = _whiten_chain_flat!(x, index, y, first(specs), hp)
    return _whiten_specs_flat!(x, index, y, Base.tail(specs), hp)
end

function _whiten_chain_flat!(x, index, y, spec::MarkovChainSpec, hp)
    free = spec.free
    if spec.centered
        n = _nfree(spec)
        # Vectorized gather (inverse of the coloring scatter): read the scattered chain values
        # `y[free.tgt]` into the contiguous latent block. Avoids scalar indexing so this also
        # raises if traced.
        @view(x[index:(index + n - 1)]) .= @view(y[free.tgt])
        return index + n
    end
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ, P = stationary_moments(p)
    @inbounds for k in 1:_nfree(spec)
        m, s = _free_moments(p, μ, P, y, free, k)
        x[index + k - 1] = (y[free.tgt[k]] - m) / s
    end
    return index + _nfree(spec)
end

function _whiten_chain_flat!(x, index, y, spec::IIDChainSpec, hp)
    @inbounds for i in spec.freeinds
        index = TV.inverse_at!(x, index, spec.fnode, y[i])
    end
    return index
end

# --- Std-space chain walkers: no Jacobian; the latent is converted to a standard
# normal via the space's cdf (exact), mirroring PT's `ScalarTransport`.

@inline _color_specs_std!(y, x, index, ::Tuple{}, hp, space) = index
@inline function _color_specs_std!(y, x, index, specs::Tuple, hp, space)
    index = _color_chain_std!(y, x, index, first(specs), hp, space)
    return _color_specs_std!(y, x, index, Base.tail(specs), hp, space)
end

function _color_chain_std!(y, x, index, spec::MarkovChainSpec, hp, space)
    # centered chains are rejected at node construction (`_check_std_transportable`)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ, P = stationary_moments(p)
    free = spec.free
    @inbounds for k in 1:_nfree(spec)
        m, s = _free_moments(p, μ, P, y, free, k)
        u = PT._clamp_unit(PT.space_cdf(space, x[index]))
        y[free.tgt[k]] = m + s * PT.space_quantile(PT.StdNormal(), u)
        index += 1
    end
    return index
end

function _color_chain_std!(y, x, index, spec::IIDChainSpec, hp, space)
    tin = PT.transport_node(spec.dist, space)
    @inbounds for i in spec.freeinds
        yi, index = PT.pfwd_step(tin, x, index)
        y[i] = yi
    end
    return index
end

@inline _whiten_specs_std!(x, index, y, ::Tuple{}, hp, space) = index
@inline function _whiten_specs_std!(x, index, y, specs::Tuple, hp, space)
    index = _whiten_chain_std!(x, index, y, first(specs), hp, space)
    return _whiten_specs_std!(x, index, y, Base.tail(specs), hp, space)
end

function _whiten_chain_std!(x, index, y, spec::MarkovChainSpec, hp, space)
    p = materialize(spec.process, _selecthp(hp, spec.hpsel))
    μ, P = stationary_moments(p)
    free = spec.free
    @inbounds for k in 1:_nfree(spec)
        m, s = _free_moments(p, μ, P, y, free, k)
        u = PT._clamp_unit(PT.space_cdf(PT.StdNormal(), (y[free.tgt[k]] - m) / s))
        x[index] = PT.space_quantile(space, u)
        index += 1
    end
    return index
end

function _whiten_chain_std!(x, index, y, spec::IIDChainSpec, hp, space)
    tin = PT.transport_node(spec.dist, space)
    @inbounds for i in spec.freeinds
        index = PT.pback_step!(x, index, tin, y[i])
    end
    return index
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
    ℓ, index = _color_specs_flat!(flag, y, x, index, values(chainspecs(d)), (;))
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
    index = _color_specs_std!(y, x, index, values(chainspecs(d)), (;), t.space)
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
    getparams(x)

Extract the instrument parameter values from a single instrument-parameter sample. For
hierarchical samples `(params = SiteArray, hyperparams = NamedTuple)` (produced by
[`GaussMarkovSitePrior`](@ref) with fitted hyperparameters) this returns `x.params`;
for plain `SiteArray` samples it is the identity.
"""
@inline getparams(x) = x
@inline getparams(x::NamedTuple{(:params, :hyperparams)}) = x.params

"""
    ObservedHierarchicalArrayPrior

The observed prior produced by an `ArrayPrior` whose [`GaussMarkovSitePrior`](@ref)
processes have fitted hyperparameters. Samples are NamedTuples
`(params = SiteArray, hyperparams = NamedTuple)` and the log-density is the exact chain
density of the parameters given the hyperparameters plus the hyperprior. Internal type.
"""
struct ObservedHierarchicalArrayPrior{D <: GaussMarkovChainDist, H <: NamedDist, S <: SiteLookup} <: Dists.ContinuousMultivariateDistribution
    dists::D
    hyperprior::H
    sitemap::S
end

Base.length(d::ObservedHierarchicalArrayPrior) = length(d.dists) + length(d.hyperprior)
Base.eltype(::ObservedHierarchicalArrayPrior) = Float64

function Dists.logpdf(d::ObservedHierarchicalArrayPrior, x::NamedTuple)
    return _chain_logpdf(d.dists, parent(x.params), x.hyperparams) +
        Dists.logpdf(d.hyperprior, x.hyperparams)
end

function Dists.rand(rng::Random.AbstractRNG, d::ObservedHierarchicalArrayPrior)
    hp = rand(rng, d.hyperprior)
    x = Vector{Float64}(undef, length(d.dists))
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
    ℓc, index = _color_specs_flat!(flag, y, x, index, values(chainspecs(d)), hp)
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
    index = _color_specs_std!(y, x, index, values(chainspecs(d)), hp, t.space)
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
    return MarkovChainSpec(sp.process, Val(hpk), inds, ts, fixedpos, sp.centered)
end

function _chainspec(sp::IIDSitePrior, site, inds, ts, fixedpos)
    free = isempty(fixedpos) ? copy(inds) : inds[setdiff(eachindex(inds), fixedpos)]
    return IIDChainSpec(sp.dist, free)
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

    finds, vals = reference_indices(array, smap, d.refant)
    fixedinds = collect(Int, finds)
    fixedvals = vals === nothing ? Float64[] : collect(float.(vals))

    # Anchored chains are exactly conditioned to zero at each site's first time, on top of
    # any reference-fixed indices (see the `anchored` docs of `GaussMarkovSitePrior`).
    for inds in lookup(smap)
        sp = getproperty(site_dists, smap.sites[first(inds)])
        (sp isa GaussMarkovSitePrior && sp.anchored) || continue
        i1 = first(inds)
        i1 ∈ fixedinds && continue
        push!(fixedinds, i1)
        push!(fixedvals, 0.0)
    end

    chains = map(lookup(smap)) do inds
        s = smap.sites[first(inds)]
        sp = getproperty(site_dists, s)
        ts = _chain_times(smap.times, inds)
        fixedpos = findall(in(fixedinds), inds)
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
