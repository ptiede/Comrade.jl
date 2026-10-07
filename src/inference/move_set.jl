# Running moves between the chunks of a sampler: the `between_chunks` hook.

export MoveSet, move_summary, move_seconds

# `counts[phase]` and `bystep[label]` are `[proposed, accepted]`; `complog` holds the log
# step scales of a componentwise move's components.
mutable struct _MoveStats
    logscale::Float64
    complog::Vector{Float64}
    counts::Dict{Symbol, Vector{Int}}
    bystep::Dict{Any, Vector{Int}}
end

function _MoveStats(m::AbstractMove)
    k = step_kind(m)
    complog = k isa ComponentwiseRandomWalk ? fill(log(k.initial_scale), k.n) : Float64[]
    return _MoveStats(
        log(_initial_scale(m)), complog, Dict(:warmup => [0, 0], :sampling => [0, 0]),
        Dict{Any, Vector{Int}}()
    )
end

"""
    step_label(move, step)

A key under which [`MoveSet`](@ref) counts the proposals and acceptances of `step`, or
`nothing` (the default) to keep only the move's totals.
"""
step_label(::AbstractMove, step) = nothing

"""
    MoveSet(post::VLBIPosterior, moves; space = nothing, θ0 = nothing, rounds = 1,
            target_accept = 0.45, output = nothing)

Metropolis–Hastings `moves` (a collection of [`AbstractMove`](@ref)s acting in the base
latent space `space` of `post`), callable as the `between_chunks` hook of a sampler:
`ms(state, tpost, info, rng) -> state`. Each call makes `rounds[j]` proposals of move `j`
(an integer applies to every move). Moves with `traceable(move)` (see
[`AbstractMove`](@ref)) go first, interleaved round by round from steps and uniforms drawn
up front; the other moves follow, interleaved the same way.

A call maps the sampler's position `z` to the base point `x = b + A z` (`A`, `b` from the
preconditioner `tpost` samples through, or the identity), runs the moves on `x`, and returns
`z = A⁻¹(x − b)`. A proposal `x′` is accepted with probability `min(1, α)`,
`log α = ℓ(x′) − ℓ(x) + log|det ∂x′/∂x|`, where `ℓ` is the log density of the base space;
the affine map has a constant Jacobian, so this is the ratio of the densities `tpost`
samples. A NaN `log α` is a rejection.

In the StdNormal space `ℓ(x) = loglik(T(x)) + log N(x; 0, I)`, so a move with `is_invariant`
is accepted on `ℓ(x′) − ℓ(x) = −(‖x′‖² − ‖x‖²)/2` without evaluating the likelihood. After
the moves, a call evaluates `ℓ` once and errors if it differs from the value the accepted
moves imply, that is, if one of them changed the likelihood. A
[`ComponentwiseRandomWalk`](@ref) move requires this acceptance: each round draws a step
for every component and accepts each component on its own, from the change of `‖x‖²` over
its `component_coords` and its log-determinant, with one step scale per component. On a device position the traceable moves' proposals,
their acceptance and the log density run as one compiled Reactant program per call (a
traced loop over rounds), and `x` stays on the device; the other moves propose on the host
from a copy of `x`. Programs are compiled once per `tpost`, so in-place refits of the
preconditioner's device buffers are seen without a recompile.

During warmup the log step scale of each [`RandomWalk`](@ref) move follows a Robbins–Monro
recursion toward `target_accept[j]`, updated after each proposal (for traceable moves, after
the call, in proposal order); it is frozen for sampling. Every warmup call logs one
line with its wall time, acceptance rates and scales; [`move_seconds`](@ref) is the total
wall time of all calls. With `output` set, the statistics of
[`move_summary`](@ref) are serialized to that path after every call.

Construction checks that the StdNormal log density has that form when any move is accepted
without the likelihood. With `θ0` (a constrained parameter point) it also checks that every
move with `is_invariant` leaves the likelihood unchanged there and that the components of
every componentwise move are independent there ([`check_components`](@ref)).
"""
struct MoveSet{M <: Tuple, C}
    moves::M
    ctx::C
    rounds::Vector{Int}
    target_accept::Vector{Float64}
    stats::Vector{_MoveStats}
    free::Vector{Bool}
    output::Union{Nothing, String}
    compiled::Base.RefValue{Any}
    seconds::Base.RefValue{Float64}
end

_per_move(v::Real, n, name) = fill(v, n)
function _per_move(v::AbstractVector, n, name)
    length(v) == n || throw(ArgumentError("$name has $(length(v)) entries for $n moves"))
    return collect(v)
end

function MoveSet(
        post::VLBIPosterior, moves; space = nothing, θ0 = nothing, rounds = 1,
        target_accept = 0.45, output = nothing
    )
    moves = Tuple(moves)
    isempty(moves) && throw(ArgumentError("MoveSet needs at least one move"))
    all(m -> m isa AbstractMove, moves) ||
        throw(ArgumentError("every move must be an AbstractMove"))
    for m in moves
        (traceable(m) && !(step_kind(m) isa RandomWalk)) &&
            throw(ArgumentError("move $(move_name(m)) is traceable but not a RandomWalk move"))
    end
    names = map(move_name, moves)
    allunique(names) || throw(ArgumentError("move names must be distinct, got $(collect(names))"))
    n = length(moves)
    r = _per_move(rounds, n, "rounds")
    all(>=(0), r) && any(>(0), r) ||
        throw(ArgumentError("rounds must be non-negative with at least one positive, got $r"))
    ta = _per_move(target_accept, n, "target_accept")
    all(a -> 0 < a < 1, ta) ||
        throw(ArgumentError("target_accept must lie in (0, 1), got $ta"))
    view = CoordinateView(post, space)
    ctx = move_context(view, moves)
    free = [_likelihood_free(view, m) for m in moves]
    for (m, f) in zip(moves, free)
        (step_kind(m) isa ComponentwiseRandomWalk && !f) && throw(
            ArgumentError(
                "move $(move_name(m)) is componentwise, which needs the StdNormal space and a " *
                    "move that leaves the likelihood unchanged"
            )
        )
    end
    any(free) && check_std_density(post, view)
    if !isnothing(θ0)
        x0 = inverse(view.tbase, θ0)
        foreach(m -> is_invariant(m) && check_invariance(m, post, view, ctx, x0), moves)
        for m in moves
            step_kind(m) isa ComponentwiseRandomWalk || continue
            check_components(m, x0, draw_step(m, Random.Xoshiro(1), _initial_scale(m)), ctx)
        end
    end
    return MoveSet(
        moves, ctx, Int.(r), Float64.(ta), [_MoveStats(m) for m in moves], free,
        isnothing(output) ? nothing : String(output), Ref{Any}(nothing), Ref(0.0)
    )
end

_likelihood_free(view::CoordinateView, m::AbstractMove) = space(view) isa PT.StdNormal && is_invariant(m)

# The change of `log N(x; 0, I)` from `x` to `x′`.
_prior_delta(x, x′) = -sum((x′ .- x) .* (x′ .+ x)) / 2

"""
    check_std_density(post, view; rtol = 1e-8)

Check that the log density of the StdNormal base space of `view` is the log-likelihood of
`post` plus `log N(x; 0, I)` up to a constant, at two points; an error otherwise.
"""
function check_std_density(post, view::CoordinateView; rtol::Real = 1.0e-8)
    tb = view.tbase
    rng = Random.Xoshiro(2)
    d(x) = logdensityof(tb, x) - loglikelihood(post, transform(tb, x)) + sum(abs2, x) / 2
    d1, d2 = d(randn(rng, dimension(tb)) ./ 2), d(randn(rng, dimension(tb)) ./ 2)
    abs(d1 - d2) <= rtol * (abs(d1) + 1) || error(
        "the StdNormal log density is not the log-likelihood plus log N(x; 0, I): their " *
            "difference changes by $(abs(d1 - d2)) between two points, so invariant moves " *
            "cannot be accepted without the likelihood"
    )
    return nothing
end

"""
    check_invariance(move, post, view, ctx, x; step, rtol = 1e-8) -> Float64

The change of the log-likelihood of `post` when `move` proposes from the base point `x` with
`step` (by default one drawn at the move's initial scale); an error if it exceeds `rtol`
relative to the log-likelihood.
"""
function check_invariance(
        move::AbstractMove, post, view::CoordinateView, ctx, x;
        step = draw_step(move, Random.Xoshiro(1), _initial_scale(move)), rtol::Real = 1.0e-8
    )
    x′, _ = propose(move, x, step, ctx)
    l0 = loglikelihood(post, transform(view.tbase, x))
    l1 = loglikelihood(post, transform(view.tbase, x′))
    abs(l1 - l0) <= rtol * (abs(l0) + 1) || error(
        "move $(move_name(move)) changed the log-likelihood from $l0 to $l1; the model is " *
            "not invariant under it"
    )
    return abs(l1 - l0)
end

# The functions a `MoveSet` call runs through: `load(position) -> z`, `to_base(z) -> x`,
# `from_base(x) -> z`, `logdensity(x) -> Float64` (of the base space),
# `prior_delta(x, x′) -> Float64`, `propose(j, x, step) -> (x′, logdet)`,
# `fused(x, ℓ, order, steps, logu) -> (x, ℓ, logα)` for the traceable moves,
# `download(x)` and `upload(xh, x)` (a host copy of `x`, and a host vector back to where `x`
# lives) and `store(z, position)`, here for host positions. The Reactant extension provides
# the device version for device positions.
function _move_kernels(ms::MoveSet, tpost, position)
    hpre() = (p = _transport_pre(tpost); isnothing(p) ? nothing : _hostify(p))
    tbase = CoordinateView(tpost.lpost, space(ms.ctx.view)).tbase
    k = (;
        load = p -> collect(Float64, vec(Array(p))),
        to_base = z -> (p = hpre(); isnothing(p) ? z : _affine_fwd(p, z)),
        from_base = x -> (p = hpre(); isnothing(p) ? x : _affine_inv(p, x)),
        logdensity = x -> Float64(logdensityof(tbase, x)),
        prior_delta = _prior_delta,
        propose = (j, x, step) -> propose(ms.moves[j], x, step, ms.ctx),
        download = copy,
        upload = (xh, x) -> xh,
        store = (z, position) -> _position_like(z, position),
    )
    tj = _traced_moves(ms)
    return merge(k, (; fused = (args...) -> _fused_sequential(k, tj, ms.free[tj], args...)))
end
_position_like(z, position) = reshape(z, size(position))

_traced_moves(ms::MoveSet) = [j for j in eachindex(ms.moves) if traceable(ms.moves[j])]

# The traceable moves' steps of one call: step `s` is a proposal of move `tj[order[s]]` with
# step `steps[s]`, accepted when `logu[s] < log α`; `free[c]` when move `tj[c]` is accepted
# without the likelihood.
function _fused_sequential(k, tj, free, x, ℓ, order, steps, logu)
    logα = similar(steps)
    for (s, c) in pairs(order)
        x′, logdet = k.propose(tj[c], x, steps[s])
        ℓ′ = free[c] ? ℓ + k.prior_delta(x, x′) : k.logdensity(x′)
        logα[s] = ℓ′ - ℓ + logdet
        if logu[s] < logα[s]
            x, ℓ = x′, ℓ′
        end
    end
    return x, ℓ, logα
end

function _kernels_for(ms::MoveSet, tpost, position)
    c = ms.compiled[]
    (!isnothing(c) && c.tpost === tpost && c.T === typeof(position)) && return c.kernels
    kernels = _move_kernels(ms, tpost, position)
    ms.compiled[] = (; tpost, T = typeof(position), kernels)
    return kernels
end

_count!(c, accepted) = (c[1] += 1; c[2] += accepted; c)

function _record_components!(st::_MoveStats, phase::Symbol, α, accepted, target, labels)
    c = st.counts[phase]
    c[1] += length(α)
    c[2] += count(accepted)
    for i in eachindex(α)
        _count!(get!(() -> [0, 0], st.bystep, labels[i]), accepted[i])
    end
    if phase === :warmup
        sweeps = st.counts[:warmup][1] ÷ length(α)
        st.complog .= clamp.(st.complog .+ (α .- target) ./ sweeps^0.6, log(1.0e-8), log(10.0))
    end
    return st
end

# One round of the componentwise move `j` on the host: a step for every component, each
# accepted on its own without the likelihood. Returns the point, its log density and the
# number of components accepted.
function _componentwise_round(ms::MoveSet, j, k, x, ℓ, phase, rng)
    m = ms.moves[j]
    st = ms.stats[j]
    xh = k.download(x)
    x′, logdet = propose(m, xh, draw_step(m, rng, exp.(st.complog)), ms.ctx)
    cs = component_coords(m)
    xn = copy(xh)
    α = zeros(length(cs))
    accepted = falses(length(cs))
    for (i, c) in pairs(cs)
        Δ = _prior_delta(xh[c], x′[c])
        a = Δ + logdet[i]
        accepted[i] = !isnan(a) && log(rand(rng)) < a
        α[i] = isnan(a) ? 0.0 : min(1.0, exp(a))
        if accepted[i]
            xn[c] = x′[c]
            ℓ += Δ
        end
    end
    _record_components!(st, phase, α, accepted, ms.target_accept[j], [component_label(m, i) for i in eachindex(cs)])
    return any(accepted) ? k.upload(xn, x) : x, ℓ, count(accepted)
end

function _record!(st::_MoveStats, kind, phase::Symbol, α, accepted::Bool, target, label)
    _count!(st.counts[phase], accepted)
    isnothing(label) || _count!(get!(() -> [0, 0], st.bystep, label), accepted)
    if phase === :warmup && kind isa RandomWalk
        n = st.counts[:warmup][1]
        st.logscale = clamp(st.logscale + (α - target) / n^0.6, log(1.0e-8), log(10.0))
    end
    return st
end

function (ms::MoveSet)(state, tpost, info, rng)
    want = space(ms.ctx.view)
    have = _base_space(tpost)
    isnothing(have) == isnothing(want) || throw(
        ArgumentError(
            "the moves act in the $(_space_name(want)) space but the sampler samples the " *
                "$(_space_name(have)) space"
        )
    )
    phase = info.phase
    phase in (:warmup, :sampling) || throw(ArgumentError("unknown sampler phase $(repr(phase))"))
    t0 = Base.time()
    position = state.position
    k = _kernels_for(ms, tpost, position)
    x = k.to_base(k.load(position))
    ℓ = k.logdensity(x)
    isfinite(ℓ) || error("the log density at the sampler's position is $ℓ before the moves")
    nacc = zeros(Int, length(ms.moves))
    tj = _traced_moves(ms)
    if !isempty(tj)
        order, steps, logu = Int[], Float64[], Float64[]
        for r in 1:maximum(ms.rounds[tj]), (c, j) in pairs(tj)
            r <= ms.rounds[j] || continue
            push!(order, c)
            push!(steps, draw_step(ms.moves[j], rng, exp(ms.stats[j].logscale)))
            push!(logu, log(rand(rng)))
        end
        x, ℓ, logα = k.fused(x, ℓ, order, steps, logu)
        for (s, c) in pairs(order)
            j = tj[c]
            m = ms.moves[j]
            a = logα[s]
            accepted = !isnan(a) && logu[s] < a
            α = isnan(a) ? 0.0 : min(1.0, exp(a))
            _record!(ms.stats[j], step_kind(m), phase, α, accepted, ms.target_accept[j], step_label(m, steps[s]))
            nacc[j] += accepted
        end
    end
    for r in 1:maximum(ms.rounds), (j, m) in enumerate(ms.moves)
        (r <= ms.rounds[j] && !(j in tj)) || continue
        st = ms.stats[j]
        kind = step_kind(m)
        if kind isa ComponentwiseRandomWalk
            x, ℓ, a = _componentwise_round(ms, j, k, x, ℓ, phase, rng)
            nacc[j] += a
            continue
        end
        step = draw_step(m, rng, exp(st.logscale))
        x′, logdet = k.propose(j, x, step)
        ℓ′ = ms.free[j] ? ℓ + k.prior_delta(x, x′) : k.logdensity(x′)
        logα = ℓ′ - ℓ + logdet
        accepted = !isnan(logα) && log(rand(rng)) < logα
        α = isnan(logα) ? 0.0 : min(1.0, exp(logα))
        _record!(st, kind, phase, α, accepted, ms.target_accept[j], step_label(m, step))
        if accepted
            x, ℓ = x′, ℓ′
            nacc[j] += 1
        end
    end
    if any(j -> ms.free[j] && nacc[j] > 0, eachindex(ms.moves))
        ℓfull = k.logdensity(x)
        # Tight: a likelihood change under this tolerance is accepted as if it were none.
        abs(ℓfull - ℓ) <= 1.0e-9 * (abs(ℓfull) + 1) || error(
            "the log density after the moves is $ℓfull, but accepting them without the " *
                "likelihood gave $ℓ: one of the accepted moves " *
                "($(join([move_name(m) for (j, m) in enumerate(ms.moves) if ms.free[j] && nacc[j] > 0], ", "))) " *
                "changed the likelihood"
        )
    end
    isnothing(ms.output) || serialize(ms.output, move_summary(ms))
    any(>(0), nacc) && (state.position = k.store(k.from_base(x), position))
    dt = Base.time() - t0
    ms.seconds[] += dt
    if phase === :warmup
        @info "moves after warmup step $(info.step)/$(info.total) ($(@sprintf("%.2f", dt)) s)\n" *
            _moves_table(ms, nacc)
    end
    return state
end

# One row per move: proposals accepted in this call, proposals made, the acceptance rate,
# and the random-walk step scale (the geometric mean over components for a componentwise
# move; blank for discrete moves).
function _moves_table(ms::MoveSet, nacc)
    rows = map(enumerate(ms.moves)) do (j, m)
        kind = step_kind(m)
        st = ms.stats[j]
        nprop = kind isa ComponentwiseRandomWalk ? ms.rounds[j] * kind.n : ms.rounds[j]
        τ = kind isa ComponentwiseRandomWalk ? exp(sum(st.complog) / length(st.complog)) :
            kind isa RandomWalk ? exp(st.logscale) : nothing
        return [
            move_name(m), string(nacc[j]), string(nprop),
            @sprintf("%.0f%%", 100 * nacc[j] / nprop),
            isnothing(τ) ? "" : @sprintf("%.3g", τ),
        ]
    end
    return pretty_table(
        String, permutedims(reduce(hcat, rows));
        column_labels = ["move", "accepted", "proposed", "rate", "step scale"],
        alignment = [:l, :r, :r, :r, :r]
    )
end

"""
    move_seconds(ms::MoveSet) -> Float64

The total wall time of the calls of `ms` so far, in seconds.
"""
move_seconds(ms::MoveSet) = ms.seconds[]

"""
    move_summary(ms::MoveSet) -> Vector{NamedTuple}

Per move: `name`, the current step scale `τ` (`nothing` for discrete moves, a vector over
the components of a componentwise move), `rounds`, and the proposals and acceptances in
warmup and sampling (a componentwise round proposes every component), plus `bystep`, the
counts over both phases per `step_label` (per `component_label` for a componentwise move)
as `label => [proposed, accepted]`.
"""
function move_summary(ms::MoveSet)
    return map(enumerate(ms.moves)) do (j, m)
        st = ms.stats[j]
        (;
            name = move_name(m),
            τ = _summary_scale(step_kind(m), st),
            rounds = ms.rounds[j],
            warmup = (proposed = st.counts[:warmup][1], accepted = st.counts[:warmup][2]),
            sampling = (proposed = st.counts[:sampling][1], accepted = st.counts[:sampling][2]),
            bystep = copy(st.bystep),
        )
    end
end

_summary_scale(::RandomWalk, st) = exp(st.logscale)
_summary_scale(::ComponentwiseRandomWalk, st) = exp.(st.complog)
_summary_scale(kind, st) = nothing
