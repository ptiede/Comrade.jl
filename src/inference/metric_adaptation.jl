# How a sampler adapts its metric during warmup. A diagonal mass matrix and a low-rank
# affine reparameterization of the latent space are two answers to the same question —
# what shape is the posterior — so they are two implementations of one interface rather
# than a flag plus a caller-supplied refit hook. A sampler asks the adaptor for its refit
# schedule, feeds it (draw, score) pairs from warmup, and asks for a new transform
# at each scheduled step.

export WelfordDiagonal, FixedMetric, FisherLowRank, GaussNewtonLowRank, NutpieDiagonal, StiffCore

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
    StiffCore(; draws = 30, directions = 8, steps = 40, power = 16, level = 1, from = 0.3)

A correction made on top of each [`FisherLowRank`](@ref) fit from the fraction `from` of
warmup on, for stiff directions that change from draw to draw. A fit to averaged draws and
scores spreads such a direction over many weak ones, so the step size stays at its stability
limit wherever the chain meets it.

At up to `draws` of the fit's window draws, a `steps`-step Lanczos run on the Hessian of the
potential in the fitted coordinates gives its `directions` largest Ritz pairs. Each draw's
Ritz pairs, written on the span `W` of all the Ritz vectors, stand for its Hessian there, and
the power mean of those matrices, `(mean_i Hᵢ^power)^(1/power)`, is diagonalized. A large
`power` follows the stiffest draw in each direction. A second Lanczos run at each draw, on the
Hessian restricted to the complement of `W`, gives the largest curvature outside it, and the
same power mean of those values is the bulk curvature `κ`. Every eigendirection with
`ψ > level κ` is whitened until `ψ = level κ`, so the core no longer sets the step size.
`level < 1` leaves margin for draws stiffer than the envelope, but narrows the core at every
draw, which slows the chain through it where it is soft. Each correction costs
`2 draws × steps` Hessian–vector products (two gradients each).
"""
struct StiffCore
    draws::Int
    directions::Int
    steps::Int
    power::Float64
    level::Float64
    from::Float64
    function StiffCore(draws, directions, steps, power, level, from)
        draws >= 2 || throw(ArgumentError("draws must be at least 2, got $draws"))
        1 <= directions <= steps ||
            throw(ArgumentError("directions must lie in 1:steps = 1:$steps, got $directions"))
        power >= 1 || throw(ArgumentError("power must be at least 1, got $power"))
        level > 0 || throw(ArgumentError("level must be positive, got $level"))
        0 <= from < 1 || throw(ArgumentError("from must lie in [0, 1), got $from"))
        return new(draws, directions, steps, power, level, from)
    end
end
StiffCore(;
    draws::Integer = 30, directions::Integer = 8, steps::Integer = 40, power::Real = 16,
    level::Real = 1, from::Real = 0.3
) = StiffCore(draws, directions, steps, power, level, from)

"""
    FisherLowRank(; rank = 16, schedule = :stan, cutoff = 2.0, γ = 1e-5,
                  min_draws = 12, max_fit_draws = 192, discard = 0.5, carry = :none,
                  max_wide = nothing, holdout = 0.0, rank_tol = 0.05, record_every = 2,
                  stiff = nothing)

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
  - `max_wide`: the most wide directions a fit keeps; `nothing` sets no limit and `0` corrects
    only stiff ones. A wide direction is one the window's draws spread far along, which a
    chain that moves slowly along it measures poorly; stretching it where the posterior is
    narrower than the window suggested makes a new stiff direction.
  - `holdout`: when positive, each fit chooses its own rank, at most `rank`, from the
    window's last `holdout` fraction of draws: a transform fit to the rest is truncated to
    the shortest prefix of its directions whose leapfrog stiffness on the held-out scores
    is within `rank_tol` of the best, and the whole window is refit at that rank. Requires
    `carry = :none` and no `fixed` directions.
  - `rank_tol`: with `holdout`, the smallest rank whose predicted leapfrog step is within
    this fraction of the best rank's is chosen.
  - `record_every`: warmup steps between recorded draws. Each NUTS step costs a whole
    trajectory of gradients while recording one costs a host copy, so frequent recording
    is nearly free. At most `2 ⌈max_fit_draws / 0.75⌉` draws are stored: past that the
    store keeps every other draw and records half as often, so a long warmup stays evenly
    covered in bounded memory (see [`FisherAdaptation`](@ref)).
  - `γ`: weight of the sample outer product against the identity in the projected draw and
    score covariances, formed as `XXᵀ/γ + I`. A smaller γ trusts the window more. Shrinking
    toward the identity rather than toward zero is what leaves a direction the window
    carries no information about at eigenvalue 1, outside the two-sided filter, instead of
    admitting it as an extreme stiff or wide direction.
  - `min_draws`: refits below this many recorded draws are skipped.
  - `max_fit_draws`: draws are thinned to at most this many columns per fit.
  - `discard`: fraction in `[0, 1)` of the draws recorded so far that each fit drops from the
    start of warmup, so fits use only the most recent draws. The first 3 draws are always
    dropped, and at least `min_draws` are always kept. The default of 0.5 keeps about the
    draws since the previous refit of a doubling schedule, as Stan's windowed adaptation
    does.
  - `carry`: what each fit builds on. `:none` fits from scratch. `:rescale` keeps the
    diagonal and directions of the transform the chain is sampling through, re-scales the
    directions on the current draws and scores, and spends the rest of `rank` on new
    directions outside their span, so directions accumulate across fits. `:keep` keeps the
    transform warmup started in (a `transport_method` preconditioner) exactly — its
    directions, scales and diagonal — and fits up to `rank` minus its rank new directions
    outside its span, replacing those of the previous fit. Use `:keep` for a starting
    transform the window cannot measure as well, such as one built from the likelihood's
    curvature.
  - `fixed`: `nothing`, or a pair `(W, sd)` of base-flat directions (the columns of `W`)
    and the standard deviation along each that every fit gives the transform. Each fit
    re-expresses them in its own diagonal standardization, takes their span with the
    variance `sd` implies along each column, and fits up to `rank` minus their number new
    directions outside that span. Use it for a direction a window cannot measure, such as
    a slow mode whose extent is known from elsewhere. Requires `carry = :none`.
  - `stiff`: `nothing`, or a [`StiffCore`](@ref) correction made on top of each fit from
    `stiff.from` of warmup on.
    The sampler supplies the Hessian–vector products (the Reactant NUTS sampler does).

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
    fixed::Union{Nothing, Tuple{Matrix{Float64}, Vector{Float64}}}
    max_wide::Union{Nothing, Int}
    holdout::Float64
    rank_tol::Float64
    record_every::Int
    stiff::Union{Nothing, StiffCore}
    function FisherLowRank{S}(
            rank, schedule, cutoff, γ, min_draws, max_fit_draws, discard, carry, fixed = nothing,
            max_wide = nothing, holdout = 0.0, rank_tol = 0.05, record_every = 2, stiff = nothing
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
        fixed = _check_fixed(fixed, rank, carry)
        (isnothing(max_wide) || max_wide >= 0) ||
            throw(ArgumentError("max_wide must be nonnegative or nothing, got $max_wide"))
        0 <= holdout < 1 || throw(ArgumentError("holdout must lie in [0, 1), got $holdout"))
        0 < rank_tol < 1 || throw(ArgumentError("rank_tol must lie in (0, 1), got $rank_tol"))
        holdout > 0 && (carry !== :none || !isnothing(fixed)) && throw(
            ArgumentError("holdout rank selection needs carry = :none and no fixed directions")
        )
        record_every >= 1 || throw(ArgumentError("record_every must be positive, got $record_every"))
        return new{S}(
            rank, schedule, cutoff, γ, min_draws, max_fit_draws, discard, carry, fixed, max_wide,
            holdout, rank_tol, record_every, stiff
        )
    end
end

function FisherLowRank(;
        rank::Int = 16, schedule = :stan, cutoff::Real = 2.0, γ::Real = 1.0e-5,
        min_draws::Int = 12, max_fit_draws::Int = 192, discard::Real = 0.5,
        carry::Symbol = :none, fixed = nothing, max_wide = nothing, holdout::Real = 0.0,
        rank_tol::Real = 0.05, record_every::Int = 2, stiff = nothing
    )
    return FisherLowRank{typeof(schedule)}(
        rank, schedule, cutoff, γ, min_draws, max_fit_draws, discard, carry, fixed, max_wide,
        holdout, rank_tol, record_every, stiff
    )
end

_check_fixed(::Nothing, rank, carry) = nothing
function _check_fixed((W, sd), rank, carry)
    W = Matrix{Float64}(W)
    sd = Vector{Float64}(sd)
    carry === :none ||
        throw(ArgumentError("fixed directions need carry = :none, got carry = $(repr(carry))"))
    size(W, 2) == length(sd) || throw(
        DimensionMismatch("fixed has $(size(W, 2)) directions but $(length(sd)) standard deviations")
    )
    0 < length(sd) < rank ||
        throw(ArgumentError("fixed needs between 1 and rank - 1 = $(rank - 1) directions, got $(length(sd))"))
    all(>(0), sd) || throw(ArgumentError("fixed standard deviations must be positive, got $sd"))
    all(c -> norm(c) > 0, eachcol(W)) || throw(ArgumentError("fixed directions must be nonzero"))
    return (W, sd)
end

function Base.show(io::IO, a::FisherLowRank)
    fx = isnothing(a.fixed) ? "" : ", $(length(a.fixed[2])) fixed directions"
    so = isnothing(a.max_wide) ? "" : a.max_wide == 0 ? ", stiff only" : ", at most $(a.max_wide) wide"
    ho = a.holdout > 0 ? ", rank from a $(a.holdout) holdout" : ""
    st = isnothing(a.stiff) ? "" : ", stiff core"
    return print(io, "FisherLowRank(rank = $(a.rank), schedule = $(repr(a.schedule))$fx$so$ho$st)")
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

Warmup draws and scores stored for [`FisherLowRank`](@ref), both in base-flat coordinates so
they stay comparable across refits, with the record number of each. Every `stride`-th
recorded draw is stored. Draws older than the next fit's window (see the adaptor's
`discard`) are dropped as records arrive, and once more than the adaptor's capacity remain,
every other stored draw is dropped and `stride` doubles. The store thus always spans the
window evenly.

[`save_adaptation!`](@ref) checkpoints it and [`load_adaptation`](@ref) restores it.
"""
mutable struct FisherAdaptation
    const draws::Vector{Vector{Float64}}
    const scores::Vector{Vector{Float64}}
    const record::Vector{Int}
    nrecords::Int
    stride::Int
end
FisherAdaptation() = FisherAdaptation(Vector{Float64}[], Vector{Float64}[], Int[], 0, 1)
FisherAdaptation(draws, scores) =
    FisherAdaptation(draws, scores, collect(eachindex(draws)), length(draws), 1)

Base.show(io::IO, st::FisherAdaptation) =
    print(io, "FisherAdaptation($(length(st.draws)) of $(st.nrecords) recorded draws, stride $(st.stride))")

_store_capacity(a::FisherLowRank) = 2 * ceil(Int, a.max_fit_draws / 0.75)

# The first record number the next fit can use after `N` records: the leading `a.discard`
# fraction (and at least the first 3 records) dropped. The cut never moves backward as
# records arrive, so draws before it are never needed again.
_window_start(a::FisherLowRank, N::Int) =
    max(round(Int, min(0.2, 3 / max(N, 1)) * N), floor(Int, a.discard * N)) + 1

# Store the draw of record `r` (already on the stride), then drop what the window no longer
# covers, keeping at least `a.min_draws`, and halve the store while it is over capacity.
function _push_draw!(a::FisherLowRank, st::FisherAdaptation, r::Int, x, g)
    push!(st.draws, x)
    push!(st.scores, g)
    push!(st.record, r)
    return _prune!(a, st)
end

function _prune!(a::FisherLowRank, st::FisherAdaptation)
    cut = _window_start(a, st.nrecords)
    ndrop = min(count(<(cut), st.record), length(st.record) - a.min_draws)
    if ndrop > 0
        deleteat!(st.draws, 1:ndrop)
        deleteat!(st.scores, 1:ndrop)
        deleteat!(st.record, 1:ndrop)
    end
    while length(st.draws) > _store_capacity(a)
        keep = reverse(lastindex(st.draws):-2:firstindex(st.draws))
        for v in (st.draws, st.scores, st.record)
            keepat!(v, keep)
        end
        st.stride *= 2
    end
    return st
end

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

function observe_draw!(a::FisherLowRank, st::FisherAdaptation, pre, position, gradient)
    x = _baseflat_draw(pre, vec(position))
    st.nrecords += 1
    if st.nrecords % st.stride == 0
        _push_draw!(a, st, st.nrecords, x, _baseflat_score(pre, vec(gradient)))
    else
        _prune!(a, st)
    end
    return x
end

"""
    adaptation_checkpoint(path, adaptor)

The checkpoint [`save_adaptation!`](@ref) writes `adaptor`'s state to. For most adaptors it
is `path` itself, rewritten in full at each save. For [`FisherLowRank`](@ref) it is a
`FisherCheckpoint`: a snapshot of the store at `path` and, at `path * ".log"`, the draws
stored since, appended at each save so that saving costs only the new draws.
"""
adaptation_checkpoint(path::AbstractString, ::AbstractMetricAdaptor) = String(path)
adaptation_checkpoint(path::AbstractString, ::FisherLowRank) = FisherCheckpoint(String(path), -1, 0, 0)

mutable struct FisherCheckpoint
    const path::String
    last::Int
    stride::Int
    logged::Int
end

_adaptation_log(path::AbstractString) = path * ".log"

"""
    save_adaptation!(checkpoint, state)

Persist the adaptation `state` to `checkpoint` (see [`adaptation_checkpoint`](@ref)).
"""
save_adaptation!(path::AbstractString, st) = (serialize(path, st); nothing)

# A new snapshot (emptying the log) is written when there is none yet, when the store has
# thinned since the last one, or when the log would hold more than half as many draws as
# the store; otherwise the draws stored since the last save are appended. The snapshot is
# replaced atomically, and the loader skips logged draws the snapshot already holds, so an
# interrupted save never loses or duplicates a draw.
function save_adaptation!(ck::FisherCheckpoint, st::FisherAdaptation)
    new = findall(>(ck.last), st.record)
    isempty(new) && ck.last >= 0 && return nothing
    log = _adaptation_log(ck.path)
    if ck.last < 0 || st.stride != ck.stride || 2 * (ck.logged + length(new)) > length(st.draws)
        tmp = ck.path * ".tmp"
        serialize(tmp, st)
        mv(tmp, ck.path; force = true)
        n = isempty(st.draws) ? 0 : length(first(st.draws))
        open(io -> write(io, Int64(n)), log, "w")
        ck.logged = 0
        ck.stride = st.stride
    else
        open(log, "a") do io
            for i in new
                write(io, Int64(st.record[i]), st.draws[i], st.scores[i])
            end
        end
        ck.logged += length(new)
    end
    ck.last = isempty(st.record) ? 0 : last(st.record)
    return nothing
end

"""
    load_adaptation(path, adaptor)

The adaptation state checkpointed at `path` by [`save_adaptation!`](@ref). For
[`FisherLowRank`](@ref) the logged draws are replayed onto the snapshot through the same
pruning and thinning the live store applied.
"""
load_adaptation(path::AbstractString, ::AbstractMetricAdaptor) = deserialize(path)

function load_adaptation(path::AbstractString, a::FisherLowRank)
    st = deserialize(path)::FisherAdaptation
    log = _adaptation_log(path)
    isfile(log) || return st
    newest = isempty(st.record) ? 0 : last(st.record)
    open(log) do io
        n = Int(read(io, Int64))
        isempty(st.draws) || length(first(st.draws)) == n || throw(
            DimensionMismatch("the adaptation log at $log holds $n-dimensional draws, the snapshot $(length(first(st.draws)))-dimensional ones")
        )
        bytes = 8 + 16 * n
        while !eof(io)
            if filesize(log) - position(io) < bytes
                @warn "Ignoring a partially written draw at the end of $log"
                break
            end
            r = Int(read(io, Int64))
            x = read!(io, Vector{Float64}(undef, n))
            g = read!(io, Vector{Float64}(undef, n))
            r <= newest && continue
            st.nrecords = r
            _push_draw!(a, st, r, x, g)
        end
    end
    return st
end

"""
    resume_adaptation!(adaptor, state, warmup_done)

Bring a loaded adaptation `state` up to a warmup resumed after `warmup_done` steps.
"""
resume_adaptation!(::AbstractMetricAdaptor, st, warmup_done::Int) = st

function resume_adaptation!(a::FisherLowRank, st::FisherAdaptation, warmup_done::Int)
    st.nrecords = warmup_done ÷ a.record_every
    return _prune!(a, st)
end

"""
    restart_step_size(adaptor, last::Bool) -> Bool

Whether step-size dual averaging restarts at a refit, `last` saying whether it is warmup's
final one. Under the `:nutpie` schedule it carries across refits, as nuts-rs adapts the step
size continuously while the metric changes, and restarts only at the last refit, for the
step-size-only end of warmup. Otherwise every refit restarts it, as Stan does at the end of
each window.
"""
restart_step_size(::AbstractMetricAdaptor, last::Bool) = true
restart_step_size(a::FisherLowRank, last::Bool) = a.schedule !== :nutpie || last

"""
    fixed_rank_cap(adaptor) -> Int or nothing

The device block's rank capacity when it is sized once for the whole warmup, or `nothing`
to grow it with the fits. The `:nutpie` schedule refits too often for each rank change to
pay a recompile, so its block holds the adaptor's full `rank`.
"""
fixed_rank_cap(::AbstractMetricAdaptor) = nothing
fixed_rank_cap(a::FisherLowRank) = a.schedule === :nutpie ? max(a.rank, 16) : nothing

"""
    record_interval(adaptor) -> Int or nothing

Warmup steps between the draws the sampler records for `adaptor`, or `nothing` to record
once per warmup chunk.
"""
record_interval(::AbstractMetricAdaptor) = nothing
record_interval(a::FisherLowRank) = a.record_every

# Positions of the `n` stored draws a fit uses, thinned evenly to at most `a.max_fit_draws`.
_refit_selection(a::FisherLowRank, n::Int) =
    unique(round.(Int, range(1, n; length = min(a.max_fit_draws, n))))

"""
    metric_refit(adaptor, state; current = nothing, initial = nothing, step = nothing,
                 n_adapts = nothing, chunk = nothing) -> LowRankPreconditioner or nothing

Fit a new latent space from the recorded draws and scores, or `nothing` when too few have
been recorded to fit from. `current` and `initial` are the preconditioners the chain samples
through now and started warmup in (`nothing` for the base space), which an adaptor with
`carry = :rescale` or `:keep` builds on. `step` is the warmup step of the refit, one of
[`metric_refit_steps`](@ref)`(adaptor, n_adapts, chunk)`, for an adaptor whose refits differ
by step.
"""
metric_refit(::AbstractMetricAdaptor, state; kwargs...) = nothing

# `p` without the zero-padded columns (`s = 1`) of a device rank slot.
function _active_directions(p::LowRankPreconditioner)
    h = _hostify(p)
    act = findall(!=(1), h.s)
    return LowRankPreconditioner(h.b, h.d, h.V[:, act], h.s[act])
end

function metric_refit(a::FisherLowRank, st::FisherAdaptation; current = nothing, initial = nothing, kwargs...)
    N = length(st.draws)
    N >= a.min_draws || return nothing
    sel = _refit_selection(a, N)
    Z = reduce(hcat, @view st.draws[sel])
    G = reduce(hcat, @view st.scores[sel])
    from = a.carry === :rescale ? current : a.carry === :keep ? initial : nothing
    carry = isnothing(from) ? nothing : _active_directions(from)
    rank = a.holdout > 0 ? _holdout_rank(a, st, 1:N) : a.rank
    return _fisher_lowrank(
        Z, G; rank, cutoff = a.cutoff, γ = a.γ, carry, keep_carried = a.carry === :keep,
        fixed = a.fixed, max_wide = a.max_wide
    )
end

"""
    _holdout_rank(a::FisherLowRank, st::FisherAdaptation, window::AbstractUnitRange) -> Int

The rank the recorded draws `window` support. A transform is fit, as `a` fits, to the first
`1 - a.holdout` of the window, and the remaining draws measure the leapfrog stiffness of each
truncation of it to its first `r` directions: the largest eigenvalue of the second moment of
their scores in the truncated transform's latent space. Scores alone measure it, so draws
that have not yet spread along a slow direction do not bias it. The stable step scales as
that eigenvalue^(-1/2); the result is the smallest `r` whose step is within `a.rank_tol` of
the best truncation's, or `a.rank` when the window is too short to split.
"""
function _holdout_rank(a::FisherLowRank, st::FisherAdaptation, window::AbstractUnitRange)
    nheld = round(Int, a.holdout * length(window))
    nfit = length(window) - nheld
    if nheld < 8 || nfit < a.min_draws
        @info "held-out rank selection skipped: $(length(window)) draws, too few to split"
        return a.rank
    end
    fitpart = first(window):(first(window) + nfit - 1)
    fsel = unique(round.(Int, range(first(fitpart), last(fitpart); length = min(a.max_fit_draws, nfit))))
    q = _fisher_lowrank(
        reduce(hcat, @view st.draws[fsel]), reduce(hcat, @view st.scores[fsel]);
        rank = a.rank, cutoff = a.cutoff, γ = a.γ, max_wide = a.max_wide
    )
    isempty(q.s) && return a.rank
    A = reduce(hcat, @view st.scores[(last(fitpart) + 1):last(window)]) .* q.d
    P = q.V' * A
    # The latent scores of the first r directions are A + Σₖ (sₖ - 1) vₖ Pₖ, with orthonormal
    # vₖ, so their Gram matrix gains (sₖ² - 1) PₖᵀPₖ per direction.
    K = A' * A
    top = [eigmax(Symmetric(K))]
    for k in eachindex(q.s)
        K .+= (q.s[k]^2 - 1) .* (P[k, :] * P[k, :]')
        push!(top, eigmax(Symmetric(K)))
    end
    best = minimum(top)
    r = findfirst(t -> sqrt(best / t) >= 1 - a.rank_tol, top) - 1
    @info(
        "held-out rank selection: $r of $(length(q.s)) directions",
        step_vs_diagonal = round(sqrt(top[1] / top[r + 1]); digits = 2),
        best_vs_diagonal = round(sqrt(top[1] / best); digits = 2),
        fit_draws = length(fsel), held_out = nheld
    )
    return r
end

"""
    GaussNewtonLowRank(curvature; rank, oversample = 10, probes_per_draw = rank + oversample,
                       threshold = 100.0, schedule = :stan, min_draws = 4, seed = 1,
                       rows = nothing, basis = nothing, residual_at = Float64[],
                       residual_draws = 30, residual_steps = 25, residual_directions = 10,
                       residual_threshold = 4.0, rescale = nothing)

Adapt the latent space to the likelihood's curvature: at each recorded warmup draw `x`,
apply the Gauss–Newton (or any symmetric positive semidefinite) curvature `H(x)` of the
negative log-likelihood to `probes_per_draw` columns of a fixed Gaussian probe matrix `Ω`
(`n × (rank + oversample)`, entries `N(0, 1/n)`, cycling through its columns), and average
the products. At each scheduled refit the averaged sketch `Ȳ ≈ H̄Ω` gives the Nyström
eigendecomposition `H̄ ≈ V Λ Vᵀ`, and the eigenpairs with `λ ≥ threshold` set a
[`LowRankPreconditioner`](@ref) with unit diagonal and scales `s = 1/√(1 + λ)`: the
posterior Hessian `I + H̄` whitened on `span(V)`. The sketch is then reset, so each fit uses
only the draws since the previous one.

Works in the `StdNormal` latent space only, where the prior is exactly N(0, I) and `I + H`
is the posterior's Gauss–Newton Hessian.

# Arguments
  - `curvature(x, W) -> Matrix`: `H(x) * W` for a base-space draw `x` (a `Vector`) and a
    matrix of directions `W`.
  - `rank`: the most eigenpairs a fit keeps. A fit with more than `rank` eigenvalues
    `≥ threshold` keeps the `rank` largest and logs how many it dropped and the largest
    dropped eigenvalue; the dropped directions stay unwhitened.
  - `oversample`: extra probe columns beyond `rank`.
  - `probes_per_draw`: probe columns applied per recorded draw, in `1:(rank + oversample)`.
    Fewer than all of them means each column averages over a different subset of draws.
  - `schedule`, `min_draws`: as for [`FisherLowRank`](@ref); a refit is also skipped until
    every probe column has been applied at least once since the last fit.
  - `seed`: the random seed of `Ω`.
  - `rows`: the latent coordinates the directions may use (strictly increasing), or
    `nothing` for all. With `rows`, `Ω` is zero off `rows`, the sketch holds only those
    rows of `H̄Ω`, and a fit's `V` is a [`RowSupportedMatrix`](@ref) on them: the
    eigendecomposition is that of `H̄` restricted to `rows`, and applying the
    preconditioner costs `length(rows)` rather than `n` per direction.
  - `basis(x) -> Matrix`: a probe matrix `Ω` built from a draw `x` in place of the Gaussian
    one, for a curvature whose range is known to lie (nearly) in a subspace that depends on
    the point. Each window between refits takes its first recorded draw as its anchor and
    probes with `basis(anchor)` at every recorded draw, applying all of its columns (so
    `rank`, `oversample` and `probes_per_draw` only bound the fit, not the probes); the
    refit is the Nyström fit on those columns. Exclusive with `rows`.
  - `residual_at`: fractions of warmup at which a residual refit runs in place of a Nyström
    fit; they must not coincide with the steps of `schedule`. A residual refit keeps the
    preconditioner the chain samples through (which must have zero shift and unit
    diagonal, as a Gauss–Newton fit has) and whitens the curvature it leaves at the draws
    recorded since the previous refit. At up to `residual_draws` of those draws, spread
    evenly, a `residual_steps`-step Lanczos run on the curvature in the current
    coordinates, `A (I + H(x)) A`, gives its `residual_directions` largest Ritz vectors.
    The curvature averaged over those draws on the span of all of them is diagonalized, and
    its directions with averaged eigenvalue `≥ residual_threshold` are whitened on top of
    the current preconditioner. The result is the symmetric square root of the combined
    covariance, so it can hold up to `residual_draws × residual_directions` directions
    beyond the current rank. This catches stiffness that no averaged sketch resolves: the
    largest eigenvectors of `H(x)` turning slightly from draw to draw, which leaves stiff
    directions outside the fitted span at draws the fit did not see. The draws must come
    after the fit: at the draws a fit was made from, it leaves little curvature. A residual
    refit costs `residual_draws × (residual_steps + r)` curvature products, `r` the number
    of collected directions.
  - `rescale`: `nothing`, or a [`NutpieDiagonal`](@ref) whose estimate sets the scale of
    each fitted direction `v` at every refit from the window's draws `x` and base-space
    gradients `g`, `σ̂ = (Var vᵀx / Var vᵀg)^(1/4)`, shrunk toward the curvature scale `s`
    and bounded to within a factor `bound` of it as that type does per coordinate. The
    curvature gives the posterior's width at a draw; where the posterior curves, the draws
    spread further along a direction than that width, and these scales follow the spread.
    It uses the same kept draws as a residual refit (and their gradients).

The sketch accumulator holds an `n × (rank + oversample)` `Float64` matrix and is
checkpointed with the warmup state, as are the window's draws kept for a residual refit
(between `residual_draws` and `2 residual_draws` vectors of length `n`). `Ω` is never
stored: each column is generated from `seed` and its index when it is used, and a draw's
products are formed 64 columns at a time. A refit works in the accumulator's memory,
generating `Ω` for the refit only and allocating the kept eigenvectors, and a failed refit
leaves the accumulator overwritten.
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
    rows::Union{Nothing, Vector{Int}}
    basis::Any
    residual_at::Vector{Float64}
    residual_draws::Int
    residual_steps::Int
    residual_directions::Int
    residual_threshold::Float64
    rescale::Any
    function GaussNewtonLowRank{F, S}(
            curvature, rank, oversample, probes_per_draw, threshold, schedule, min_draws, seed,
            rows = nothing, basis = nothing, residual_at = Float64[], residual_draws = 30, residual_steps = 25,
            residual_directions = 10, residual_threshold = 4.0, rescale = nothing
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
        isnothing(rows) || (!isempty(rows) && issorted(rows; lt = <=) && first(rows) >= 1) ||
            throw(ArgumentError("rows must be a non-empty, strictly increasing list of coordinates"))
        (isnothing(rows) || isnothing(basis)) ||
            throw(ArgumentError("rows and basis are exclusive: a basis sets its own coordinates"))
        all(f -> 0 < f < 1, residual_at) ||
            throw(ArgumentError("residual_at must hold fractions of warmup in (0, 1), got $residual_at"))
        residual_draws >= 1 || throw(ArgumentError("residual_draws must be positive, got $residual_draws"))
        1 <= residual_directions <= residual_steps || throw(
            ArgumentError(
                "residual_directions must lie in 1:residual_steps = 1:$residual_steps, got $residual_directions"
            )
        )
        residual_threshold > 1 ||
            throw(ArgumentError("residual_threshold must exceed 1 (unit curvature is already white), got $residual_threshold"))
        (isnothing(rescale) || rescale isa NutpieDiagonal) ||
            throw(ArgumentError("rescale must be nothing or a NutpieDiagonal, got $(typeof(rescale))"))
        (isnothing(rescale) || min_draws >= 2) ||
            throw(ArgumentError("rescale estimates variances and needs min_draws ≥ 2, got $min_draws"))
        return new{F, S}(
            curvature, rank, oversample, probes_per_draw, threshold, schedule, min_draws, seed,
            isnothing(rows) ? nothing : collect(Int, rows), basis, sort(collect(Float64, residual_at)),
            residual_draws, residual_steps,
            residual_directions, residual_threshold, rescale
        )
    end
end

function GaussNewtonLowRank(
        curvature; rank::Integer, oversample::Integer = 10,
        probes_per_draw::Integer = rank + oversample, threshold::Real = 100.0,
        schedule = :stan, min_draws::Integer = 4, seed::Integer = 1, rows = nothing,
        basis = nothing, residual_at = Float64[], residual_draws::Integer = 30, residual_steps::Integer = 25,
        residual_directions::Integer = 10, residual_threshold::Real = 4.0, rescale = nothing
    )
    return GaussNewtonLowRank{typeof(curvature), typeof(schedule)}(
        curvature, rank, oversample, probes_per_draw, threshold, schedule, min_draws, seed, rows,
        basis, residual_at, residual_draws, residual_steps, residual_directions, residual_threshold,
        rescale
    )
end

function Base.show(io::IO, a::GaussNewtonLowRank)
    res = (isempty(a.residual_at) ? "" : ", residual_at = $(a.residual_at)") *
        (isnothing(a.rescale) ? "" : ", rescale = $(a.rescale)")
    return print(
        io, "GaussNewtonLowRank(rank = $(a.rank), threshold = $(a.threshold), " *
            "schedule = $(repr(a.schedule))$res)"
    )
end

check_metric_space(::GaussNewtonLowRank, ::PT.StdNormal) = nothing
check_metric_space(::GaussNewtonLowRank, space) = throw(
    ArgumentError(
        "GaussNewtonLowRank needs the StdNormal latent space (transport_method = " *
            "StdNormal() or a Preconditioned StdNormal transform), got $(_space_name(space))"
    )
)

function metric_refit_steps(a::GaussNewtonLowRank, n_adapts::Int, chunk::Int)
    fits = filter(s -> 0 < s < n_adapts, _refit_steps(a.schedule, n_adapts, chunk))
    res = _residual_steps(a, n_adapts, chunk)
    shared = intersect(fits, res)
    isempty(shared) || throw(
        ArgumentError("residual_at and schedule share the warmup steps $shared of $n_adapts")
    )
    return sort!(vcat(fits, res))
end

_residual_steps(a::GaussNewtonLowRank, n_adapts::Int, chunk::Int) =
    isempty(a.residual_at) ? Int[] :
    filter(s -> 0 < s < n_adapts, _refit_steps(a.residual_at, n_adapts, chunk))

# The coordinates the directions use, for latent dimension `n`.
function _probe_rows(a::GaussNewtonLowRank, n::Int)
    isnothing(a.rows) && return 1:n
    last(a.rows) <= n || throw(
        ArgumentError("rows reach coordinate $(last(a.rows)) of a latent space of dimension $n")
    )
    return a.rows
end

# The number of probe columns, `rank + oversample`, checked against the number of rows they
# span among the `n` latent coordinates.
function _nprobes(a::GaussNewtonLowRank, n::Int)
    k = a.rank + a.oversample
    nr = length(_probe_rows(a, n))
    k <= nr || throw(
        ArgumentError("rank + oversample = $k exceeds the $nr latent coordinates the directions may use")
    )
    return k
end

# Columns `cols` of the probe matrix `Ω` for latent dimension `n`: column `j` is a standard
# normal vector scaled by `1/√n`, generated from the adaptor's seed and `j`.
function _probe_columns(a::GaussNewtonLowRank, n::Int, cols)
    W = Matrix{Float64}(undef, n, length(cols))
    for (c, j) in pairs(cols)
        ω = view(W, :, c)
        randn!(Random.Xoshiro(hash((a.seed, j))), ω)
        ω ./= sqrt(n)
    end
    return W
end

# `Ω` in the coordinates `rows` of a latent space of dimension `n`.
_probe_matrix(a::GaussNewtonLowRank, n::Int) = _probe_columns(a, length(_probe_rows(a, n)), 1:_nprobes(a, n))

"""
    GaussNewtonSketch

The sketch accumulated by [`GaussNewtonLowRank`](@ref) since its last fit: `Y[:, j]` is the
sum of the curvature products with probe column `j` over the `counts[j]` draws that applied
it (on the adaptor's `rows`), `next` the next column to apply, `ndraws` the draws recorded
`n` the latent dimension (0 before the first draw) and `anchor` the draw the window's
probes are built from when the adaptor has a `basis` (empty otherwise, and between a fit
and the next draw). `kept` holds every `stride`-th draw of the window, and `scores` their
base-space gradients, for a residual refit or a rescale; when it reaches twice the adaptor's `residual_draws`, every other one is dropped and
`stride` doubles, so the kept draws stay evenly spread. Serializable.
"""
mutable struct GaussNewtonSketch
    Y::Matrix{Float64}
    counts::Vector{Int}
    next::Int
    ndraws::Int
    n::Int
    anchor::Vector{Float64}
    kept::Vector{Vector{Float64}}
    scores::Vector{Vector{Float64}}
    stride::Int
end
GaussNewtonSketch() = GaussNewtonSketch(
    Matrix{Float64}(undef, 0, 0), Int[], 1, 0, 0, Float64[], Vector{Float64}[], Vector{Float64}[], 1
)

Base.show(io::IO, st::GaussNewtonSketch) =
    print(io, "GaussNewtonSketch($(st.ndraws) draws, $(length(st.counts)) probes)")

init_metric_adaptation(::GaussNewtonLowRank) = GaussNewtonSketch()

function _reset!(st::GaussNewtonSketch)
    fill!(st.Y, 0)
    fill!(st.counts, 0)
    st.ndraws = 0
    empty!(st.anchor)
    empty!(st.kept)
    empty!(st.scores)
    st.stride = 1
    return st
end

# Keep `x` and its score `g`, the window's draw number `st.ndraws + 1`, when it falls on the
# stride.
function _keep_draw!(a::GaussNewtonLowRank, st::GaussNewtonSketch, x, g)
    (isempty(a.residual_at) && isnothing(a.rescale)) && return st
    if st.ndraws % st.stride == 0
        push!(st.kept, copy(x))
        push!(st.scores, copy(g))
    end
    if length(st.kept) >= 2 * a.residual_draws
        st.kept = st.kept[1:2:end]
        st.scores = st.scores[1:2:end]
        st.stride *= 2
    end
    return st
end

function observe_draw!(a::GaussNewtonLowRank, st::GaussNewtonSketch, pre, position, gradient)
    x = _baseflat_draw(pre, vec(position))
    _keep_draw!(a, st, x, _baseflat_score(pre, vec(gradient)))
    isnothing(a.basis) || return _observe_basis!(a, st, x)
    n = length(x)
    k = _nprobes(a, n)
    nr = length(_probe_rows(a, n))
    if isempty(st.counts)
        st.Y = zeros(nr, k)
        st.counts = zeros(Int, k)
        st.n = n
    end
    size(st.Y) == (nr, k) || throw(
        DimensionMismatch("the sketch holds $(size(st.Y)) products, but the draw has $nr probed coordinates and $k probes")
    )
    cols = mod1.(st.next .+ (0:(a.probes_per_draw - 1)), k)
    _apply_probes!(a, st, x, cols)
    st.counts[cols] .+= 1
    st.next = mod1(st.next + a.probes_per_draw, k)
    st.ndraws += 1
    return x
end

# A recorded draw for an adaptor with a `basis`: the window's first draw becomes its anchor,
# and every column of `basis(anchor)` is applied at `x`.
function _observe_basis!(a::GaussNewtonLowRank, st::GaussNewtonSketch, x)
    n = length(x)
    isempty(st.anchor) && append!(st.anchor, x)
    Ω = _basis_matrix(a, st.anchor)
    k = size(Ω, 2)
    if size(st.Y) != (n, k)
        st.Y = zeros(n, k)
        st.counts = zeros(Int, k)
    end
    st.n = n
    for b in Iterators.partition(1:k, 64)
        HΩ = a.curvature(x, Ω[:, b])
        size(HΩ) == (n, length(b)) || throw(
            DimensionMismatch("curvature returned a $(size(HΩ)) matrix for $(length(b)) directions of dimension $n")
        )
        all(isfinite, HΩ) || error(
            "the curvature product at a warmup draw is not finite (draw $(st.ndraws + 1) since the last fit)"
        )
        st.Y[:, b] .+= HΩ
    end
    st.counts .+= 1
    st.ndraws += 1
    return x
end

function _basis_matrix(a::GaussNewtonLowRank, anchor)
    Ω = a.basis(anchor)
    (size(Ω, 1) == length(anchor) && size(Ω, 2) >= 1) || throw(
        DimensionMismatch("basis returned a $(size(Ω)) matrix for a latent dimension of $(length(anchor))")
    )
    return Matrix{Float64}(Ω)
end

# Add the curvature products at `x` with probe columns `cols` to the sketch, `block` columns
# per call of the curvature, so only one `n × block` slice of probes and products exists at a
# time. With `rows`, the probes are embedded in the latent space and only those rows of the
# products are kept.
function _apply_probes!(a::GaussNewtonLowRank, st::GaussNewtonSketch, x, cols; block::Int = 64)
    n = length(x)
    rows = _probe_rows(a, n)
    for b in Iterators.partition(cols, block)
        Ωb = _probe_columns(a, length(rows), b)
        W = isnothing(a.rows) ? Ωb : (Wf = zeros(n, length(b)); Wf[rows, :] = Ωb; Wf)
        HΩ = a.curvature(x, W)
        size(HΩ) == (n, length(b)) || throw(
            DimensionMismatch("curvature returned a $(size(HΩ)) matrix for $(length(b)) directions of dimension $n")
        )
        all(isfinite, HΩ) || error(
            "the curvature product at a warmup draw is not finite (draw $(st.ndraws + 1) since the last fit)"
        )
        st.Y[:, b] .+= isnothing(a.rows) ? HΩ : view(HΩ, rows, :)
    end
    return st
end

# The eigenvalues, decreasing, of the Nyström approximation `Y (ΩᵀY)⁻¹ Yᵀ` of a positive
# semidefinite `H` from its sketch `Y = HΩ` (Tropp, Yurtsever, Udell & Cevher 2017, Alg. 3),
# and the eigenvectors of those at least `cut`. The approximation is unchanged by `Ω → ΩR⁻¹`,
# so with `ΩᵀΩ = RᵀR` the algorithm runs on the orthonormal `Q = ΩR⁻¹` and `HQ = YR⁻¹`; the
# shift `ν` keeps the core factorization stable for a rank-deficient `H`. `Y` is
# overwritten: it becomes the shifted `HQ + νQ` and then the thin factor `B`, whose singular
# pairs come from the `k × k` Gram matrix `BᵀB`. `Ω` is generated for the call only.
function _nystrom!(a::GaussNewtonLowRank, Y::AbstractMatrix, cut)
    n, k = size(Y)
    return _nystrom!(Y, _probe_columns(a, n, 1:k), cut)
end

function _nystrom!(Y::AbstractMatrix, Ω::AbstractMatrix, cut)
    n, k = size(Y)
    R = cholesky(Symmetric(Ω' * Ω)).U
    ν = sqrt(n) * eps(norm(Y) * opnorm(inv(R)))
    Y .+= ν .* Ω
    M = Ω' * Y
    Ω = nothing
    rdiv!(Y, R)
    C = cholesky(Symmetric(R' \ M / R); check = false)
    issuccess(C) || error(
        "the Gauss–Newton sketch is not positive definite on its probes; with " *
            "probes_per_draw below rank + oversample each probe averages different draws — " *
            "raise probes_per_draw or min_draws"
    )
    rdiv!(Y, C.U)
    E = eigen(Symmetric(Y' * Y); sortby = -)
    σ² = max.(E.values, 0)
    λ = max.(σ² .- ν, 0)
    keep = findall(>=(cut), λ)
    return λ, keep, Y * (E.vectors[:, keep] ./ sqrt.(σ²[keep])')
end

function metric_refit(
        a::GaussNewtonLowRank, st::GaussNewtonSketch; current = nothing, step = nothing,
        n_adapts = nothing, chunk = nothing, kwargs...
    )
    if !isempty(a.residual_at)
        any(isnothing, (step, n_adapts, chunk)) && throw(
            ArgumentError("a Gauss–Newton adaptor with residual_at needs the refit's step, n_adapts and chunk")
        )
        step in _residual_steps(a, n_adapts, chunk) && return _residual_refit(a, st, current)
    end
    (st.ndraws >= a.min_draws && !isempty(st.counts) && all(>(0), st.counts)) || return nothing
    n = st.n
    st.Y ./= st.counts'
    λ, keep, U = isnothing(a.basis) ? _nystrom!(a, st.Y, a.threshold) :
        _nystrom!(st.Y, _basis_matrix(a, st.anchor), a.threshold)
    # `keep` is in decreasing eigenvalue order, so truncation drops the weakest directions.
    if length(keep) > a.rank
        @info "Gauss–Newton refit: $(length(keep)) eigenvalues exceed threshold = " *
            "$(a.threshold); keeping the largest rank = $(a.rank), the largest dropped is " *
            "λ = $(round(λ[keep[a.rank + 1]]; sigdigits = 4))"
        keep = keep[1:a.rank]
        U = U[:, 1:a.rank]
    end
    s = _rescale(a, st, n, U, inv.(sqrt.(1 .+ λ[keep])))
    V = isnothing(a.rows) ? U : RowSupportedMatrix(n, a.rows, U)
    _reset!(st)
    return LowRankPreconditioner(zeros(n), ones(n), V, s)
end

# A residual refit: `current` with the curvature it leaves at the window's kept draws
# whitened on top.
function _residual_refit(a::GaussNewtonLowRank, st::GaussNewtonSketch, current)
    isempty(st.kept) && return nothing
    isnothing(current) && throw(
        ArgumentError("a residual refit corrects the current preconditioner, but the chain samples without one")
    )
    p = _active_directions(current)
    (all(iszero, p.b) && all(isone, p.d)) || throw(
        ArgumentError("a residual refit needs a current preconditioner with zero shift and unit diagonal")
    )
    n = length(p.b)
    U = if isnothing(a.rows)
        Matrix(p.V)
    else
        (p.V isa RowSupportedMatrix && p.V.rows == a.rows) || throw(
            ArgumentError("a residual refit with rows needs a current preconditioner supported on those rows")
        )
        Matrix(p.V.M)
    end
    sel = unique(round.(Int, range(1, length(st.kept), length = min(a.residual_draws, length(st.kept)))))
    U, s = _residual_stage(a, st.kept[sel], n, U, collect(Float64, p.s))
    s = _rescale(a, st, n, U, s)
    _reset!(st)
    V = isnothing(a.rows) ? U : RowSupportedMatrix(n, a.rows, U)
    return LowRankPreconditioner(zeros(n), ones(n), V, s)
end

# The scales `s` of the directions `U` (on the probed rows) re-estimated from the window's
# kept draws and scores (see `GaussNewtonLowRank`).
function _rescale(a::GaussNewtonLowRank, st::GaussNewtonSketch, n, U, s)
    isnothing(a.rescale) && return s
    m = length(st.kept)
    m >= 2 || throw(ArgumentError("rescaling needs at least two kept draws, got $m"))
    rows = _probe_rows(a, n)
    vx = vec(var(U' * reduce(hcat, (x[rows] for x in st.kept)); dims = 2))
    vg = vec(var(U' * reduce(hcat, (g[rows] for g in st.scores)); dims = 2))
    c = log(a.rescale.bound)
    w = m / (m + a.rescale.n0)
    t = map(s, vx, vg) do sk, x, g
        (iszero(x) || iszero(g)) && return sk
        sk * exp(c * tanh(w * ((log(x) - log(g)) / 4 - log(sk)) / c))
    end
    f = sort(t ./ s)
    q(p) = round(f[clamp(round(Int, p * length(f)), 1, length(f))]; sigdigits = 3)
    @info "Gauss–Newton rescale from $m draws: scale factors over $(length(f)) directions " *
        "10% $(q(0.1)), median $(q(0.5)), 90% $(q(0.9)), largest $(round(last(f); sigdigits = 3))"
    return t
end

# The curvature in the coordinates `y` of the fit `U, s` (on the probed rows), at the draw
# `x`, applied to the columns of `Y`: `A (Y + H(x) A Y)` with `A = I + U (s - 1) Uᵀ`,
# `block` columns per call of the curvature.
function _fitted_curvature(a::GaussNewtonLowRank, x, n, U, s, Y; block::Int = 64)
    rows = _probe_rows(a, n)
    A(Z) = Z .+ U * ((s .- 1) .* (U' * Z))
    AY = A(Y)
    out = similar(AY)
    for b in Iterators.partition(axes(AY, 2), block)
        W = isnothing(a.rows) ? AY[:, b] : (Wf = zeros(n, length(b)); Wf[rows, :] = AY[:, b]; Wf)
        HW = a.curvature(x, W)
        size(HW) == (n, length(b)) || throw(
            DimensionMismatch("curvature returned a $(size(HW)) matrix for $(length(b)) directions of dimension $n")
        )
        all(isfinite, HW) || error("the curvature product at a kept warmup draw is not finite")
        out[:, b] .= AY[:, b] .+ (isnothing(a.rows) ? HW : view(HW, rows, :))
    end
    return A(out)
end

# The `k` largest Ritz values and vectors of `steps` Lanczos steps (full reorthogonalization) of the
# symmetric operator `op` on vectors of length `m`, from a random start.
function _lanczos_top(op, m::Int, steps::Int, k::Int, rng)
    steps = min(steps, m)
    Q = zeros(m, steps + 1)
    α, β = zeros(steps), zeros(steps)
    q = randn(rng, m)
    Q[:, 1] = q ./ norm(q)
    nsteps = steps
    for j in 1:steps
        w = vec(op(Q[:, j:j]))
        α[j] = dot(w, view(Q, :, j))
        Qj = view(Q, :, 1:j)
        for _ in 1:2
            w .-= Qj * (Qj' * w)
        end
        β[j] = norm(w)
        if β[j] <= 1.0e-10 * abs(α[j])
            nsteps = j
            break
        end
        Q[:, j + 1] = w ./ β[j]
    end
    E = eigen(SymTridiagonal(α[1:nsteps], β[1:(nsteps - 1)]); sortby = -)
    top = 1:min(k, nsteps)
    return E.values[top], Q[:, 1:nsteps] * E.vectors[:, top]
end

# The residual stage of a residual refit (see `GaussNewtonLowRank`): the directions `U` (on the probed
# rows) and scales `s` of a preconditioner with the stiff directions left in its coordinates at the
# draws `xs` whitened on top.
function _residual_stage(a::GaussNewtonLowRank, xs, n, U, s)
    m = length(_probe_rows(a, n))
    t0 = Base.time()
    ritz = [
        _lanczos_top(
                Y -> _fitted_curvature(a, x, n, U, s, Y), m, a.residual_steps, a.residual_directions,
                Random.Xoshiro(hash((a.seed, :residual, i)))
            ) for (i, x) in pairs(xs)
    ]
    peaks = sort(map(r -> first(first(r)), ritz))
    C = reduce(hcat, last.(ritz))
    F = svd(C)
    W = F.U[:, findall(>(1.0e-3 * first(F.S)), F.S)]
    M = sum(x -> W' * _fitted_curvature(a, x, n, U, s, W), xs) ./ length(xs)
    E = eigen(Symmetric(M); sortby = -)
    keep = findall(>=(a.residual_threshold), E.values)
    @info "Gauss–Newton residual refit: largest curvature left at the draws " *
        "$(round(first(peaks); sigdigits = 3))–$(round(last(peaks); sigdigits = 3)) " *
        "(median $(round(peaks[(end + 1) ÷ 2]; sigdigits = 3))); $(size(W, 2)) directions from $(length(xs)) draws, " *
        "averaged curvature up to $(round(first(E.values); sigdigits = 4)); whitening " *
        "$(length(keep)) at or above $(a.residual_threshold) ($(round(Base.time() - t0; digits = 1)) s)"
    isempty(keep) && return U, s
    return _whiten_on_top(U, s, W * E.vectors[:, keep], inv.(sqrt.(E.values[keep])))
end

# The directions and scales of the preconditioner `A = I + U (s - 1) Uᵀ` followed, in its
# coordinates, by `B = I + Wk (t - 1) Wkᵀ` (orthonormal `Wk`): the covariance `A B² A` is the
# identity off span(U, Wk), and its symmetric square root is fitted on that span.
function _whiten_on_top(U, s, Wk, t)
    Q = Matrix(qr(hcat(U, Wk)).Q)
    AQ = Q .+ U * ((s .- 1) .* (U' * Q))
    G = AQ .+ Wk * ((t .- 1) .* (Wk' * AQ))
    K = eigen(Symmetric(G' * G))
    μ = K.values
    all(>(0), μ) || error("the combined covariance is not positive definite")
    act = findall(v -> abs(sqrt(v) - 1) > 1.0e-8, μ)
    return Q * K.vectors[:, act], sqrt.(μ[act])
end

stiff_core(::AbstractMetricAdaptor) = nothing
stiff_core(a::FisherLowRank) = a.stiff

"""
    stiff_core_due(adaptor, step, n_adapts) -> Bool

Whether the refit at warmup step `step` of `n_adapts` gets the adaptor's
[`StiffCore`](@ref) correction.
"""
stiff_core_due(a::AbstractMetricAdaptor, step::Integer, n_adapts::Integer) =
    (sc = stiff_core(a); !isnothing(sc) && step >= sc.from * n_adapts)

# Orthonormal directions `D` and values `ψ > cap` of the power mean
# `(mean_i Psᵢ^power)^(1/power)` of the positive semidefinite `Ps`. The mean of powers resolves
# only values near the largest (the rest fall below its rounding error), so the directions
# within a factor 2 of the top are taken and the rest of the space is treated again,
# with the matrices projected onto it.
function _power_envelope(Ps, power, cap)
    r = size(first(Ps), 1)
    Q = Matrix{Float64}(I, r, r)
    D, ψ = zeros(r, 0), Float64[]
    while size(Q, 2) > 0
        Pq = [Symmetric(Q' * P * Q) for P in Ps]
        top = maximum(eigmax, Pq)
        top > cap || break
        E = eigen(Symmetric(sum(P -> (P ./ top)^power, Pq) ./ length(Pq)); sortby = -)
        vals = top .* max.(E.values, 0) .^ (1 / power)
        first(vals) > cap || break
        done = findall(>=(max(first(vals) / 2, cap)), vals)
        D = hcat(D, Q * E.vectors[:, done])
        append!(ψ, vals[done])
        Q = Q * E.vectors[:, setdiff(axes(E.vectors, 2), done)]
    end
    return D, ψ
end

"""
    stiff_core_refit(a::FisherLowRank, st::FisherAdaptation, pre, hvp) -> LowRankPreconditioner

`pre`, the transform fitted by `a`'s refit, with the [`StiffCore`](@ref) correction `a.stiff`
on top, from the draws of the fit's window. `hvp(x, Y)` returns the Hessian of the potential
(negative log density) at the base-flat point `x` applied to the columns of `Y`.
"""
function stiff_core_refit(a::FisherLowRank, st::FisherAdaptation, pre, hvp)
    sc = a.stiff
    p = _active_directions(pre)
    win = _refit_selection(a, length(st.draws))
    xs = st.draws[win[unique(round.(Int, range(1, length(win), length = min(sc.draws, length(win)))))]]
    U, s = Matrix(p.V), collect(Float64, p.s)
    A(Y) = Y .+ U * ((s .- 1) .* (U' * Y))
    # The Hessian in the fitted coordinates z, with x = b + d ⊙ A z: Aᵀ D H(x) D A.
    fitted(x, Y) = A(p.d .* hvp(x, p.d .* A(Y)))
    m = length(p.b)
    t0 = Base.time()
    ritz = [
        _lanczos_top(Y -> fitted(x, Y), m, sc.steps, sc.directions, Random.Xoshiro(hash((:stiff, i))))
            for (i, x) in pairs(xs)
    ]
    F = svd(reduce(hcat, last.(ritz)))
    W = F.U[:, findall(>(1.0e-3 * first(F.S)), F.S)]
    outside(Y) = Y .- W * (W' * Y)
    tops = [
        first(first(_lanczos_top(Y -> outside(fitted(x, outside(Y))), m, sc.steps, 1, Random.Xoshiro(hash((:bulk, i))))))
            for (i, x) in pairs(xs)
    ]
    pmean(v) = mean(t -> max(t, 0)^sc.power, v)^(1 / sc.power)
    κ = pmean(tops)
    # Each draw's Hessian on W from its own Ritz pairs: its stiff part lies in their span.
    Ps = map(ritz) do (vals, vecs)
        C = W' * vecs
        Symmetric(C * Diagonal(max.(vals, 0)) * C')
    end
    cap = sc.level * κ
    D, ψ = _power_envelope(Ps, sc.power, cap)
    peak = pmean(map(r -> first(first(r)), ritz))
    @info "stiff-core refit at $(length(xs)) draws: power-mean top curvature $(round(peak; sigdigits = 3)), " *
        "outside the $(size(W, 2))-direction core $(round(κ; sigdigits = 3)); whitening $(length(ψ)) " *
        "directions above $(round(cap; sigdigits = 3)) ($(round(Base.time() - t0; digits = 1)) s)"
    isempty(ψ) && return pre
    V2, s2 = _whiten_on_top(U, s, W * D, sqrt.(cap ./ ψ))
    return LowRankPreconditioner(p.b, p.d, V2, s2)
end

"""
    NutpieDiagonal(; n0 = 5, bound = 10)

A diagonal inverse mass matrix estimated from draws and gradients, as in nutpie: the
variance of coordinate `i` is `σ̂ᵢ² = √(Var xᵢ / Var gᵢ)` over the draws `x` and log-density
gradients `g` of the current warmup segment, which is exact for a Gaussian with diagonal
covariance. With
`n` draws the estimate is shrunk toward the unit scale and bounded smoothly to
`[1/bound, bound]`:

    log σᵢ = c tanh(w log σ̂ᵢ / c),   w = n / (n + n0),   c = log(bound),

and the inverse mass matrix is `σ²`. A coordinate with constant draws and gradients carries
no information and keeps the unit scale. Set as `ReactantNUTS(; diagonal)` alongside a
metric adaptor that does not adapt its own diagonal.
"""
struct NutpieDiagonal
    n0::Float64
    bound::Float64
    function NutpieDiagonal(n0, bound)
        n0 >= 0 || throw(ArgumentError("n0 must be nonnegative, got $n0"))
        bound > 1 || throw(ArgumentError("bound must exceed 1, got $bound"))
        return new(n0, bound)
    end
end
NutpieDiagonal(; n0::Real = 5, bound::Real = 10) = NutpieDiagonal(n0, bound)

# Running means and sums of squared deviations (Welford) of the draws and gradients.
mutable struct _DiagonalMoments
    n::Int
    mx::Vector{Float64}
    sx::Vector{Float64}
    mg::Vector{Float64}
    sg::Vector{Float64}
end
_DiagonalMoments(d::Integer) = _DiagonalMoments(0, zeros(d), zeros(d), zeros(d), zeros(d))

function _observe!(m::_DiagonalMoments, x::AbstractVector, g::AbstractVector)
    eachindex(x, g) == eachindex(m.mx) ||
        throw(DimensionMismatch("draw, gradient and moments must match: $(axes(x)), $(axes(g)), $(axes(m.mx))"))
    m.n += 1
    for i in eachindex(x, g, m.mx)
        δx = x[i] - m.mx[i]
        m.mx[i] += δx / m.n
        m.sx[i] += δx * (x[i] - m.mx[i])
        δg = g[i] - m.mg[i]
        m.mg[i] += δg / m.n
        m.sg[i] += δg * (g[i] - m.mg[i])
    end
    return m
end

function _reset!(m::_DiagonalMoments)
    m.n = 0
    fill!(m.mx, 0); fill!(m.sx, 0); fill!(m.mg, 0); fill!(m.sg, 0)
    return m
end

"""
    inverse_mass(d::NutpieDiagonal, m) -> Vector

The diagonal inverse mass matrix of `d` from the moments `m`; it needs at least two draws.
"""
function inverse_mass(d::NutpieDiagonal, m::_DiagonalMoments)
    m.n >= 2 || throw(ArgumentError("the diagonal estimate needs at least two draws, got $(m.n)"))
    c = log(d.bound)
    w = m.n / (m.n + d.n0)
    return map(m.sx, m.sg) do sx, sg
        (iszero(sx) && iszero(sg)) && return 1.0
        exp(2c * tanh(w * (log(sx) - log(sg)) / (4c)))
    end
end
