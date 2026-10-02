# How a sampler adapts its metric during warmup. A diagonal mass matrix and a low-rank
# affine reparameterization of the latent space are two answers to the same question —
# what shape is the posterior — so they are two implementations of one interface rather
# than a flag plus a caller-supplied refit hook. A sampler asks the adaptor for its refit
# schedule, feeds it one (draw, score) pair per warmup chunk, and asks for a new transform
# at each scheduled step.

export WelfordDiagonal, FixedMetric, FisherLowRank, GaussNewtonLowRank

"""
    AbstractMetricAdaptor

Warmup metric-adaptation strategy for a sampler. Implementations:
[`WelfordDiagonal`](@ref), [`FixedMetric`](@ref), [`FisherLowRank`](@ref),
[`GaussNewtonLowRank`](@ref).

The sampler-facing protocol is:

  - [`adapts_welford`](@ref) — does the backend's own diagonal adaptation run?
  - [`check_metric_space`](@ref) — can the adaptor work in this latent space?
  - [`metric_refit_steps`](@ref) — warmup steps at which to refit.
  - [`init_metric_adaptation`](@ref) — build the accumulator threaded through warmup.
  - [`observe_draw!`](@ref) — record one warmup draw and its score.
  - [`metric_refit`](@ref) — fit a new latent space from what has been recorded.
"""
abstract type AbstractMetricAdaptor end

"""
    WelfordDiagonal()

Adapt a diagonal mass matrix by Welford accumulation over Stan's warmup windows. This is
the sampler's built-in adaptation and the default.
"""
struct WelfordDiagonal <: AbstractMetricAdaptor end

"""
    FixedMetric()

Do not adapt the metric: it stays at the identity and only the step size adapts. Use this
when the latent space already carries the geometry — a preconditioning transform composed
into the posterior IS the metric, and diagonal adaptation would renormalize the marginals
and undo any marginal-versus-conditional balance the transform encodes.
"""
struct FixedMetric <: AbstractMetricAdaptor end

"""
    FisherLowRank(; rank = 16, schedule = :stan, cutoff = 2.0, γ = 1e-5,
                  min_draws = 12, max_fit_draws = 192, discard = 0.0, carry = :none)

Adapt the latent space itself: at each scheduled warmup step, fit a
[`LowRankPreconditioner`](@ref) by the Fisher-divergence estimator of Seyboldt, Carlson
& Carpenter (arXiv:2603.18845) from the warmup draws collected so far and their
scores, and continue warmup in the fitted coordinates. The metric stays at the identity
throughout — the transform carries the geometry.

Draws and scores come from the sampler's own state, mapped back to the base-flat space
that every transform routes through, so they describe the distribution the chain actually
explored. Fitting always happens in base-flat coordinates and replaces the transform
rather than composing onto it: low-rank-plus-diagonal maps are not closed under
composition, and base-flat is the one reference frame that is invariant across refits.

# Arguments
  - `rank`: maximum number of corrected directions per fit.
  - `schedule`: `:stan` (a first fit at step 100, then doubling gaps to 85% of warmup),
    `:nutpie` (fits from the start of warmup, every chunk to 30% and every eighth chunk
    to 85%), or a vector of fractions of warmup in `(0, 1)`.
  - `cutoff`: an eigenvalue is corrected when `λ ≥ cutoff` (wide) or `λ ≤ 1/cutoff` (stiff).
  - `γ`: weight of the sample outer product against the identity in the projected draw and
    score covariances, formed as `XXᵀ/γ + I`. A smaller γ trusts the window more. Shrinking
    toward the identity rather than toward zero is what leaves a direction the window
    carries no information about at eigenvalue 1, outside the two-sided filter, instead of
    admitting it as an extreme stiff or wide direction.
  - `min_draws`: refits below this many recorded draws are skipped.
  - `max_fit_draws`: draws are thinned to at most this many columns per fit.
  - `discard`: fraction in `[0, 1)` of the draws recorded so far that each fit drops from the
    start of warmup, so fits use only the most recent draws. The first 3 draws are always
    dropped, and at least `min_draws` are always kept.
  - `carry`: what each fit builds on. `:none` fits from scratch. `:rescale` keeps the
    directions of the transform the chain is sampling through, re-scaled on the current
    draws and scores, and spends the rest of `rank` on new directions outside their span,
    so directions accumulate across fits. `:keep` keeps the transform warmup started in
    (a `transport_method` preconditioner) exactly — its directions, scales and diagonal —
    and fits up to `rank` minus its rank new directions outside its span, replacing those
    of the previous fit. Use `:keep` for a starting transform the window cannot measure as
    well, such as one built from the likelihood's curvature.

To start warmup in fitted coordinates rather than on a unit metric, pass a
[`LowRankPreconditioner`](@ref) as the sampler's `transport_method`.
"""
struct FisherLowRank{S} <: AbstractMetricAdaptor
    rank::Int
    schedule::S
    cutoff::Float64
    γ::Float64
    min_draws::Int
    max_fit_draws::Int
    discard::Float64
    carry::Symbol
    function FisherLowRank{S}(
            rank, schedule, cutoff, γ, min_draws, max_fit_draws, discard, carry
        ) where {S}
        rank > 0 || throw(ArgumentError("rank must be positive, got $rank"))
        cutoff > 1 || throw(ArgumentError("cutoff must be greater than 1, got $cutoff"))
        γ > 0 || throw(ArgumentError("the regularization weight γ must be positive, got $γ"))
        min_draws >= 4 ||
            throw(ArgumentError("min_draws must be at least 4, got $min_draws"))
        max_fit_draws >= min_draws || throw(
            ArgumentError(
                "max_fit_draws ($max_fit_draws) must be at least min_draws ($min_draws)"
            )
        )
        0 <= discard < 1 ||
            throw(ArgumentError("discard must lie in [0, 1), got $discard"))
        _check_refit_schedule(schedule)
        carry in (:none, :rescale, :keep) ||
            throw(ArgumentError("carry must be :none, :rescale or :keep, got $(repr(carry))"))
        return new{S}(rank, schedule, cutoff, γ, min_draws, max_fit_draws, discard, carry)
    end
end

function FisherLowRank(;
        rank::Int = 16, schedule = :stan, cutoff::Real = 2.0, γ::Real = 1.0e-5,
        min_draws::Int = 12, max_fit_draws::Int = 192, discard::Real = 0.0,
        carry::Symbol = :none
    )
    return FisherLowRank{typeof(schedule)}(
        rank, schedule, cutoff, γ, min_draws, max_fit_draws, discard, carry
    )
end

function Base.show(io::IO, a::FisherLowRank)
    return print(io, "FisherLowRank(rank = $(a.rank), schedule = $(repr(a.schedule)))")
end

_check_refit_schedule(s::Symbol) = s in (:stan, :nutpie) || throw(
    ArgumentError(
        "unknown refit schedule :$s; use :stan, :nutpie, or a vector of fractions in (0, 1)"
    )
)
_check_refit_schedule(v::AbstractVector{<:Real}) = all(f -> 0 < f < 1, v) || throw(
    ArgumentError("manual refit fractions must lie strictly in (0, 1), got $v")
)
_check_refit_schedule(x) = throw(
    ArgumentError(
        "refit schedule must be :stan, :nutpie, or a vector of fractions in (0, 1), got $x"
    )
)

"""
    check_metric_space(adaptor, space)

Throw if `adaptor` cannot adapt a sampler working in front of the latent space `space`
(`nothing` for flat, or a `StdNormal`).
"""
check_metric_space(::AbstractMetricAdaptor, space) = nothing

"""
    adapts_welford(adaptor) -> Bool

Whether the sampler's own diagonal (Welford) mass-matrix adaptation should run.
"""
adapts_welford(::AbstractMetricAdaptor) = false
adapts_welford(::WelfordDiagonal) = true

"""
    metric_refit_steps(adaptor, n_adapts, chunk) -> Vector{Int}

Warmup steps at which the sampler should ask for a new latent space. Steps outside
`(0, n_adapts)` are dropped, so a short warmup simply gets fewer refits.
"""
metric_refit_steps(::AbstractMetricAdaptor, n_adapts::Int, chunk::Int) = Int[]
metric_refit_steps(a::FisherLowRank, n_adapts::Int, chunk::Int) =
    filter(s -> 0 < s < n_adapts, _refit_steps(a.schedule, n_adapts, chunk))

_refit_steps(v::AbstractVector{<:Real}, na::Int, chunk::Int) =
    sort!(unique(round.(Int, v .* na)))

function _refit_steps(s::Symbol, na::Int, chunk::Int)
    if s === :stan
        # A structural fit at step 100, then exponentially growing gaps until 85% of
        # warmup. Early fits are frequent while the geometry still moves; late segments
        # run uninterrupted for thousands of steps, which is when step-size precision
        # matters — each refit restarts dual averaging, so segment length IS the
        # stabilization budget. The gap is capped so the metric still gets late refreshes
        # from the richest windows.
        steps = Int[100]
        gap = 5 * chunk
        cap = max(round(Int, 0.3 * na), 10 * chunk)
        while last(steps) + gap <= round(Int, 0.85 * na)
            push!(steps, last(steps) + gap)
            gap = min(2 * gap, cap)
        end
        return steps
    end
    # nutpie's cadence (Seyboldt, Carlson & Carpenter 2026, §3): updates from the very
    # start of warmup — every chunk to 30% of warmup, every eighth chunk to 85%, and
    # step-size only for the tail.
    return sort!(
        unique(
            vcat(
                collect(100:chunk:round(Int, 0.3 * na)),
                collect(round(Int, 0.3 * na):(8 * chunk):round(Int, 0.85 * na)),
            )
        )
    )
end

"""
    FisherAdaptation

Warmup draws and scores accumulated for [`FisherLowRank`](@ref), both in base-flat
coordinates so they stay comparable across refits. Serializable: checkpointing it beside
the sampler state is what lets an interrupted warmup resume without losing the window it
had built up.
"""
struct FisherAdaptation
    draws::Vector{Vector{Float64}}
    scores::Vector{Vector{Float64}}
end
FisherAdaptation() = FisherAdaptation(Vector{Float64}[], Vector{Float64}[])

Base.show(io::IO, st::FisherAdaptation) =
    print(io, "FisherAdaptation($(length(st.draws)) draws)")

"""
    init_metric_adaptation(adaptor) -> state

Accumulator threaded through warmup, or `nothing` when the adaptor keeps no state.
"""
init_metric_adaptation(::AbstractMetricAdaptor) = nothing
init_metric_adaptation(::FisherLowRank) = FisherAdaptation()

# The preconditioner a transformed posterior is currently sampling through, or `nothing`
# when it samples its base space (flat or StdNormal) directly.
_transport_pre(tpost::TransformedVLBIPosterior) =
    _node_pre(PT.transport_node(tpost.transform))
_node_pre(node::PreconditionedFlat) = node.pre
_node_pre(node::PreconditionedStd) = node.pre
_node_pre(::Any) = nothing

_baseflat_draw(::Nothing, z::AbstractVector) = collect(Float64, z)
_baseflat_draw(pre, z::AbstractVector) = Float64.(_affine_fwd(_pre_for(pre, z), z))
_baseflat_score(::Nothing, g::AbstractVector) = collect(Float64, g)
_baseflat_score(pre, g::AbstractVector) = Float64.(_affine_invT(_pre_for(pre, g), g))

"""
    observe_draw!(adaptor, state, pre, position, gradient) -> Vector or nothing

Record one warmup draw in `state`, the accumulator of `adaptor`. `position` and `gradient` are the sampler's own, expressed in the
latent space the preconditioner `pre` defines (`nothing` for plain base-flat); both are
mapped back to base-flat coordinates — the draw through the transform, the gradient
through its inverse transpose — so that draws from different refit segments describe one
distribution. Returns the base-flat draw, which is also the position a new transform must
be re-expressed from.
"""
observe_draw!(::AbstractMetricAdaptor, ::Nothing, pre, position, gradient) = nothing

function observe_draw!(::FisherLowRank, st::FisherAdaptation, pre, position, gradient)
    x = _baseflat_draw(pre, vec(position))
    push!(st.draws, x)
    push!(st.scores, _baseflat_score(pre, vec(gradient)))
    return x
end

# Indices of the `N` recorded draws a fit uses: the first `a.discard` fraction (and at least the
# first 3 draws) dropped, at least `a.min_draws` kept, and the rest thinned to at most
# `a.max_fit_draws`.
function _refit_selection(a::FisherLowRank, N::Int)
    ndrop = max(round(Int, min(0.2, 3 / N) * N), floor(Int, a.discard * N))
    start = min(ndrop, N - a.min_draws) + 1
    return unique(
        round.(Int, range(start, N, length = min(a.max_fit_draws, N - start + 1)))
    )
end

"""
    metric_refit(adaptor, state; current = nothing, initial = nothing)
        -> LowRankPreconditioner or nothing

Fit a new latent space from the recorded draws and scores, or `nothing` when too few have
been recorded to fit from. `current` and `initial` are the preconditioners the chain samples
through now and started warmup in (`nothing` for the base space), which an adaptor with
`carry = :rescale` or `:keep` builds on.
"""
metric_refit(::AbstractMetricAdaptor, state; kwargs...) = nothing

# `p` without the zero-padded columns (`s = 1`) of a device rank slot.
function _active_directions(p::LowRankPreconditioner)
    h = _hostify(p)
    act = findall(!=(1), h.s)
    return LowRankPreconditioner(h.b, h.d, h.V[:, act], h.s[act])
end

function metric_refit(a::FisherLowRank, st::FisherAdaptation; current = nothing, initial = nothing)
    N = length(st.draws)
    N >= a.min_draws || return nothing
    sel = _refit_selection(a, N)
    Z = reduce(hcat, @view st.draws[sel])
    G = reduce(hcat, @view st.scores[sel])
    from = a.carry === :rescale ? current : a.carry === :keep ? initial : nothing
    carry = isnothing(from) ? nothing : _active_directions(from)
    return _fisher_lowrank(
        Z, G; rank = a.rank, cutoff = a.cutoff, γ = a.γ, carry, keep_carried = a.carry === :keep
    )
end

"""
    GaussNewtonLowRank(curvature; rank, oversample = 10, probes_per_draw = rank + oversample,
                       threshold = 100.0, schedule = :stan, min_draws = 4, seed = 1)

Adapt the latent space to the likelihood's curvature: at each recorded warmup draw `x`,
apply the Gauss–Newton (or any symmetric positive semidefinite) curvature `H(x)` of the
negative log-likelihood to `probes_per_draw` columns of a fixed random orthonormal probe
matrix `Ω` (`n × (rank + oversample)`, cycling through its columns), and average the
products. At each scheduled refit the averaged sketch `Ȳ ≈ H̄Ω` gives the Nyström
eigendecomposition `H̄ ≈ V Λ Vᵀ`, and the eigenpairs with `λ ≥ threshold` set a
[`LowRankPreconditioner`](@ref) with unit diagonal and scales `s = 1/√(1 + λ)`: the
posterior Hessian `I + H̄` whitened on `span(V)`. The sketch is then reset, so each fit uses
only the draws since the previous one.

Works in the `StdNormal` latent space only, where the prior is exactly N(0, I) and `I + H`
is the posterior's Gauss–Newton Hessian.

# Arguments
  - `curvature(x, W) -> Matrix`: `H(x) * W` for a base-space draw `x` (a `Vector`) and a
    matrix of directions `W`.
  - `rank`: the most eigenpairs a fit may keep. A fit with more than `rank` eigenvalues
    `≥ threshold` errors (the spectrum would be truncated; raise `rank`).
  - `oversample`: extra probe columns beyond `rank`.
  - `probes_per_draw`: probe columns applied per recorded draw, in `1:(rank + oversample)`.
    Fewer than all of them means each column averages over a different subset of draws.
  - `schedule`, `min_draws`: as for [`FisherLowRank`](@ref); a refit is also skipped until
    every probe column has been applied at least once since the last fit.
  - `seed`: the random seed of `Ω`.

The sketch accumulator holds an `n × (rank + oversample)` `Float64` matrix and is
checkpointed with the warmup state.
"""
struct GaussNewtonLowRank{F, S} <: AbstractMetricAdaptor
    curvature::F
    rank::Int
    oversample::Int
    probes_per_draw::Int
    threshold::Float64
    schedule::S
    min_draws::Int
    seed::Int
    probes::Base.RefValue{Matrix{Float64}}
    function GaussNewtonLowRank{F, S}(
            curvature, rank, oversample, probes_per_draw, threshold, schedule, min_draws, seed
        ) where {F, S}
        rank > 0 || throw(ArgumentError("rank must be positive, got $rank"))
        oversample >= 1 ||
            throw(ArgumentError("oversample must be at least 1 to detect a truncated spectrum, got $oversample"))
        1 <= probes_per_draw <= rank + oversample || throw(
            ArgumentError("probes_per_draw must lie in 1:$(rank + oversample) (rank + oversample), got $probes_per_draw")
        )
        threshold > 0 || throw(ArgumentError("threshold must be positive, got $threshold"))
        min_draws >= 1 || throw(ArgumentError("min_draws must be at least 1, got $min_draws"))
        _check_refit_schedule(schedule)
        return new{F, S}(
            curvature, rank, oversample, probes_per_draw, threshold, schedule, min_draws, seed,
            Ref(Matrix{Float64}(undef, 0, 0))
        )
    end
end

function GaussNewtonLowRank(
        curvature; rank::Integer, oversample::Integer = 10,
        probes_per_draw::Integer = rank + oversample, threshold::Real = 100.0,
        schedule = :stan, min_draws::Integer = 4, seed::Integer = 1
    )
    return GaussNewtonLowRank{typeof(curvature), typeof(schedule)}(
        curvature, rank, oversample, probes_per_draw, threshold, schedule, min_draws, seed
    )
end

function Base.show(io::IO, a::GaussNewtonLowRank)
    return print(
        io, "GaussNewtonLowRank(rank = $(a.rank), threshold = $(a.threshold), " *
            "schedule = $(repr(a.schedule)))"
    )
end

check_metric_space(::GaussNewtonLowRank, ::PT.StdNormal) = nothing
check_metric_space(::GaussNewtonLowRank, space) = throw(
    ArgumentError(
        "GaussNewtonLowRank needs the StdNormal latent space (transport_method = " *
            "StdNormal() or a Preconditioned StdNormal transform), got $(_space_name(space))"
    )
)

metric_refit_steps(a::GaussNewtonLowRank, n_adapts::Int, chunk::Int) =
    filter(s -> 0 < s < n_adapts, _refit_steps(a.schedule, n_adapts, chunk))

# The probe matrix `Ω` for latent dimension `n`, generated from the adaptor's seed on first
# use, so the serialized sketch need not carry it.
function _probe_matrix(a::GaussNewtonLowRank, n::Int)
    k = a.rank + a.oversample
    k <= n || throw(
        ArgumentError("rank + oversample = $k exceeds the latent dimension $n")
    )
    if size(a.probes[]) != (n, k)
        a.probes[] = _thinq(randn(Random.Xoshiro(a.seed), n, k), k)
    end
    return a.probes[]
end

"""
    GaussNewtonSketch

The sketch accumulated by [`GaussNewtonLowRank`](@ref) since its last fit: `Y[:, j]` is the
sum of the curvature products with probe column `j` over the `counts[j]` draws that applied
it, `next` the next column to apply and `ndraws` the draws recorded. Serializable.
"""
mutable struct GaussNewtonSketch
    Y::Matrix{Float64}
    counts::Vector{Int}
    next::Int
    ndraws::Int
end
GaussNewtonSketch() = GaussNewtonSketch(Matrix{Float64}(undef, 0, 0), Int[], 1, 0)

Base.show(io::IO, st::GaussNewtonSketch) =
    print(io, "GaussNewtonSketch($(st.ndraws) draws, $(length(st.counts)) probes)")

init_metric_adaptation(::GaussNewtonLowRank) = GaussNewtonSketch()

function _reset!(st::GaussNewtonSketch)
    fill!(st.Y, 0)
    fill!(st.counts, 0)
    st.ndraws = 0
    return st
end

function observe_draw!(a::GaussNewtonLowRank, st::GaussNewtonSketch, pre, position, gradient)
    x = _baseflat_draw(pre, vec(position))
    n = length(x)
    Ω = _probe_matrix(a, n)
    k = size(Ω, 2)
    if isempty(st.counts)
        st.Y = zeros(n, k)
        st.counts = zeros(Int, k)
    end
    size(st.Y) == (n, k) || throw(
        DimensionMismatch("the sketch holds $(size(st.Y)) products, but the draw has dimension $n and $k probes")
    )
    cols = mod1.(st.next .+ (0:(a.probes_per_draw - 1)), k)
    HΩ = a.curvature(x, Ω[:, cols])
    size(HΩ) == (n, length(cols)) || throw(
        DimensionMismatch("curvature returned a $(size(HΩ)) matrix for $(length(cols)) directions of dimension $n")
    )
    all(isfinite, HΩ) || error(
        "the curvature product at a warmup draw is not finite (draw $(st.ndraws + 1) since the last fit)"
    )
    st.Y[:, cols] .+= HΩ
    st.counts[cols] .+= 1
    st.next = mod1(st.next + a.probes_per_draw, k)
    st.ndraws += 1
    return x
end

# Eigenpairs of the Nyström approximation `Y (ΩᵀY)⁻¹ Yᵀ` of a positive semidefinite `H` from
# its sketch `Y = HΩ` (Tropp, Yurtsever, Udell & Cevher 2017, Alg. 3), eigenvalues
# decreasing. The shift `ν` keeps the core factorization stable for a rank-deficient `H`.
function _nystrom(Ω::AbstractMatrix, Y::AbstractMatrix)
    ν = sqrt(size(Y, 1)) * eps(norm(Y))
    Yν = Y .+ ν .* Ω
    C = cholesky(Symmetric(Ω' * Yν); check = false)
    issuccess(C) || error(
        "the Gauss–Newton sketch is not positive definite on its probes; with " *
            "probes_per_draw below rank + oversample each probe averages different draws — " *
            "raise probes_per_draw or min_draws"
    )
    F = svd(Yν / C.U)
    return max.(F.S .^ 2 .- ν, 0), F.U
end

function metric_refit(a::GaussNewtonLowRank, st::GaussNewtonSketch; kwargs...)
    (st.ndraws >= a.min_draws && !isempty(st.counts) && all(>(0), st.counts)) || return nothing
    n = size(st.Y, 1)
    λ, U = _nystrom(_probe_matrix(a, n), st.Y ./ st.counts')
    keep = findall(>=(a.threshold), λ)
    length(keep) <= a.rank || error(
        "$(length(keep)) Gauss–Newton eigenvalues exceed threshold = $(a.threshold), more " *
            "than rank = $(a.rank): the spectrum would be truncated; raise rank or threshold"
    )
    _reset!(st)
    return LowRankPreconditioner(zeros(n), ones(n), U[:, keep], inv.(sqrt.(1 .+ λ[keep])))
end
