# How a sampler adapts its metric during warmup. A diagonal mass matrix and a low-rank
# affine reparameterization of the latent space are two answers to the same question —
# what shape is the posterior — so they are two implementations of one interface rather
# than a flag plus a caller-supplied refit hook. A sampler asks the adaptor for its refit
# schedule, feeds it one (draw, score) pair per warmup chunk, and asks for a new transform
# at each scheduled step.

export WelfordDiagonal, FixedMetric, FisherLowRank

"""
    AbstractMetricAdaptor

Warmup metric-adaptation strategy for a sampler. Implementations:
[`WelfordDiagonal`](@ref), [`FixedMetric`](@ref), [`FisherLowRank`](@ref).

The sampler-facing protocol is:

  - [`adapts_welford`](@ref) — does the backend's own diagonal adaptation run?
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
                  s_floor = 0.02, min_draws = 12, max_fit_draws = 192)

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
  - `γ`: ridge added to the projected draw and score covariances. With far fewer draws
    than dimensions, sampling noise in the joint draw-plus-score subspace otherwise passes
    the two-sided filter as spurious directions.
  - `s_floor`: lower bound on a fitted direction scale. An underflowing stiff scale
    compresses its axis by `1/s` and drags the step size to zero with it.
  - `min_draws`: refits below this many recorded draws are skipped.
  - `max_fit_draws`: draws are thinned to at most this many columns per fit.

To start warmup in fitted coordinates rather than on a unit metric, pass a
[`LowRankPreconditioner`](@ref) as the sampler's `transport_method`.
"""
struct FisherLowRank{S} <: AbstractMetricAdaptor
    rank::Int
    schedule::S
    cutoff::Float64
    γ::Float64
    s_floor::Float64
    min_draws::Int
    max_fit_draws::Int
    function FisherLowRank{S}(
            rank, schedule, cutoff, γ, s_floor, min_draws, max_fit_draws
        ) where {S}
        rank > 0 || throw(ArgumentError("rank must be positive, got $rank"))
        cutoff > 1 || throw(ArgumentError("cutoff must be greater than 1, got $cutoff"))
        γ > 0 || throw(ArgumentError("the ridge γ must be positive, got $γ"))
        s_floor > 0 || throw(ArgumentError("s_floor must be positive, got $s_floor"))
        min_draws >= 4 ||
            throw(ArgumentError("min_draws must be at least 4, got $min_draws"))
        max_fit_draws >= min_draws || throw(
            ArgumentError(
                "max_fit_draws ($max_fit_draws) must be at least min_draws ($min_draws)"
            )
        )
        _check_refit_schedule(schedule)
        return new{S}(rank, schedule, cutoff, γ, s_floor, min_draws, max_fit_draws)
    end
end

function FisherLowRank(;
        rank::Int = 16, schedule = :stan, cutoff::Real = 2.0, γ::Real = 1.0e-5,
        s_floor::Real = 0.02, min_draws::Int = 12, max_fit_draws::Int = 192
    )
    return FisherLowRank{typeof(schedule)}(
        rank, schedule, cutoff, γ, s_floor, min_draws, max_fit_draws
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
# when its latent space is plain base-flat.
_transport_pre(tpost::TransformedVLBIPosterior) =
    _node_pre(PT.transport_node(tpost.transform))
_node_pre(node::PreconditionedFlat) = node.pre
_node_pre(::Any) = nothing

_baseflat_draw(::Nothing, z::AbstractVector) = collect(Float64, z)
_baseflat_draw(pre, z::AbstractVector) = Float64.(_affine_fwd(_pre_for(pre, z), z))
_baseflat_score(::Nothing, g::AbstractVector) = collect(Float64, g)
_baseflat_score(pre, g::AbstractVector) = Float64.(_affine_invT(_pre_for(pre, g), g))

"""
    observe_draw!(state, pre, position, gradient) -> Vector or nothing

Record one warmup draw. `position` and `gradient` are the sampler's own, expressed in the
latent space the preconditioner `pre` defines (`nothing` for plain base-flat); both are
mapped back to base-flat coordinates — the draw through the transform, the gradient
through its inverse transpose — so that draws from different refit segments describe one
distribution. Returns the base-flat draw, which is also the position a new transform must
be re-expressed from.
"""
observe_draw!(::Nothing, pre, position, gradient) = nothing

function observe_draw!(st::FisherAdaptation, pre, position, gradient)
    x = _baseflat_draw(pre, vec(position))
    push!(st.draws, x)
    push!(st.scores, _baseflat_score(pre, vec(gradient)))
    return x
end

"""
    metric_refit(adaptor, state) -> LowRankPreconditioner or nothing

Fit a new latent space from the recorded draws and scores, or `nothing` when too few have
been recorded to fit from.
"""
metric_refit(::AbstractMetricAdaptor, state) = nothing

function metric_refit(a::FisherLowRank, st::FisherAdaptation)
    N = length(st.draws)
    N >= a.min_draws || return nothing
    # Drop the true transient only — a run started from a near-posterior draw has a short
    # one, and a fractional discard would starve the early fit windows.
    start = round(Int, min(0.2, 3 / N) * N) + 1
    sel = unique(
        round.(Int, range(start, N, length = min(a.max_fit_draws, N - start + 1)))
    )
    Z = reduce(hcat, @view st.draws[sel])
    G = reduce(hcat, @view st.scores[sel])
    return _fisher_lowrank(
        Z, G; rank = a.rank, cutoff = a.cutoff, γ = a.γ, s_floor = a.s_floor
    )
end
