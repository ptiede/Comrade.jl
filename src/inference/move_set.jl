# Running moves between the chunks of a sampler: the `between_chunks` hook.

export MoveSet, move_summary

mutable struct _MoveStats
    logscale::Float64
    proposed::Dict{Symbol, Int}
    accepted::Dict{Symbol, Int}
    bystep::Dict{Any, Vector{Int}}
end

_MoveStats(m::AbstractMove) = _MoveStats(
    log(_initial_scale(step_kind(m))), Dict(:warmup => 0, :sampling => 0),
    Dict(:warmup => 0, :sampling => 0), Dict{Any, Vector{Int}}()
)
_initial_scale(k::RandomWalk) = k.initial_scale
_initial_scale(::DiscreteSymmetric) = 1.0

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
(an integer applies to every move), interleaved round by round.

A proposal maps the base point `x = b + A z` of the sampler's position `z` (`A`, `b` from
the preconditioner `tpost` samples through, or the identity) to `x′`, sets
`z′ = z + A⁻¹(x′ − x)`, and accepts with probability `min(1, α)`,
`log α = ℓ(z′) − ℓ(z) + log|det ∂x′/∂x|`, where `ℓ` is the log density of `tpost`; a NaN
`log α` is a rejection. On a device position `ℓ` is a compiled Reactant program (compiled
once per `tpost`); proposals are made on the host.

During warmup the log step scale of each [`RandomWalk`](@ref) move follows a Robbins–Monro
recursion toward `target_accept[j]`; it is frozen for sampling. Every warmup call logs one
line of acceptance rates and scales. With `output` set, the statistics of
[`move_summary`](@ref) are serialized to that path after every call.

With `θ0` (a constrained parameter point), construction checks that every move with
`is_invariant` leaves the likelihood unchanged there.
"""
struct MoveSet{M <: Tuple, C}
    moves::M
    ctx::C
    rounds::Vector{Int}
    target_accept::Vector{Float64}
    stats::Vector{_MoveStats}
    output::Union{Nothing, String}
    compiled::Base.RefValue{Any}
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
    if !isnothing(θ0)
        x0 = inverse(view.tbase, θ0)
        foreach(m -> is_invariant(m) && check_invariance(m, post, view, ctx, x0), moves)
    end
    return MoveSet(
        moves, ctx, Int.(r), Float64.(ta), [_MoveStats(m) for m in moves],
        isnothing(output) ? nothing : String(output), Ref{Any}(nothing)
    )
end

"""
    check_invariance(move, post, view, ctx, x; rtol = 1e-8)

Error unless a proposal of `move` from the base point `x` leaves the log-likelihood of
`post` unchanged (up to `rtol`).
"""
function check_invariance(move::AbstractMove, post, view::CoordinateView, ctx, x; rtol::Real = 1.0e-8)
    kind = step_kind(move)
    x′, _ = propose(move, x, draw_step(move, Random.Xoshiro(1), _initial_scale(kind)), ctx)
    l0 = loglikelihood(post, transform(view.tbase, x))
    l1 = loglikelihood(post, transform(view.tbase, x′))
    abs(l1 - l0) <= rtol * (abs(l0) + 1) || error(
        "move $(move_name(move)) changed the log-likelihood from $l0 to $l1; the model is " *
            "not invariant under it"
    )
    return nothing
end

# A host function of a host latent vector giving the log density of `tpost`. The Reactant
# extension compiles it when the sampler's position lives on the device.
_logdensity_closure(tpost, position) = z -> logdensityof(tpost, z)
_position_like(z, position) = reshape(z, size(position))

function _logdensity_for(ms::MoveSet, tpost, position)
    c = ms.compiled[]
    (!isnothing(c) && c.tpost === tpost) && return c.f
    f = _logdensity_closure(tpost, position)
    ms.compiled[] = (; tpost, f)
    return f
end

_base_point(::Nothing, z) = z
_base_point(pre, z) = _affine_fwd(_hostify(pre), z)

# The sampled-space step that moves the base point by `Δx`: the linear part of the inverse
# preconditioner, so coordinates a move does not touch stay exactly where they were when
# no preconditioner is composed.
_sampled_step(::Nothing, Δx) = Δx
function _sampled_step(pre, Δx)
    p = _hostify(pre)
    w = Δx ./ p.d
    isempty(p.s) && return w
    return w .+ p.V * ((inv.(p.s) .- 1) .* (p.V' * w))
end

function _record!(st::_MoveStats, kind, phase::Symbol, α, accepted::Bool, target, label)
    st.proposed[phase] += 1
    st.accepted[phase] += accepted
    if !isnothing(label)
        c = get!(() -> [0, 0], st.bystep, (phase, label))
        c[1] += 1
        c[2] += accepted
    end
    if phase === :warmup && kind isa RandomWalk
        n = st.proposed[:warmup]
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
    position = state.position
    z = collect(Float64, vec(Array(position)))
    pre = _transport_pre(tpost)
    ℓf = _logdensity_for(ms, tpost, position)
    ℓ = ℓf(z)
    isfinite(ℓ) || error("the log density at the sampler's position is $ℓ before the moves")
    x = _base_point(pre, z)
    nacc = zeros(Int, length(ms.moves))
    for r in 1:maximum(ms.rounds), (j, m) in enumerate(ms.moves)
        r <= ms.rounds[j] || continue
        st = ms.stats[j]
        kind = step_kind(m)
        step = draw_step(m, rng, exp(st.logscale))
        x′, logdet = propose(m, x, step, ms.ctx)
        z′ = z .+ _sampled_step(pre, x′ .- x)
        ℓ′ = ℓf(z′)
        logα = ℓ′ - ℓ + logdet
        accepted = !isnan(logα) && log(rand(rng)) < logα
        α = isnan(logα) ? 0.0 : min(1.0, exp(logα))
        _record!(st, kind, phase, α, accepted, ms.target_accept[j], step_label(m, step))
        if accepted
            x, z, ℓ = x′, z′, ℓ′
            nacc[j] += 1
        end
    end
    if phase === :warmup
        parts = map(enumerate(ms.moves)) do (j, m)
            s = "$(move_name(m)) $(nacc[j])/$(ms.rounds[j])"
            step_kind(m) isa RandomWalk ? s * " τ=$(@sprintf("%.3g", exp(ms.stats[j].logscale)))" : s
        end
        @info "moves after warmup step $(info.step)/$(info.total): " * join(parts, "; ")
    end
    isnothing(ms.output) || serialize(ms.output, move_summary(ms))
    state.position = _position_like(z, position)
    return state
end

"""
    move_summary(ms::MoveSet) -> Vector{NamedTuple}

Per move: `name`, the current step scale `τ` (`nothing` for discrete moves), `rounds`, and
the proposals and acceptances in warmup and sampling, plus `bystep`, the counts per
`step_label` as `(phase, label) => [proposed, accepted]`.
"""
function move_summary(ms::MoveSet)
    return map(enumerate(ms.moves)) do (j, m)
        st = ms.stats[j]
        (;
            name = move_name(m),
            τ = step_kind(m) isa RandomWalk ? exp(st.logscale) : nothing,
            rounds = ms.rounds[j],
            warmup = (proposed = st.proposed[:warmup], accepted = st.accepted[:warmup]),
            sampling = (proposed = st.proposed[:sampling], accepted = st.accepted[:sampling]),
            bystep = copy(st.bystep),
        )
    end
end
