# Running moves on a sampler's position: `MoveSet`, its statistics and its joint-distribution
# test.

export MoveSet, move_summary, move_seconds, check_moveset

# `logscales` holds the log step scales: none for a move without one, one per component for a
# componentwise move. `counts[phase]` and `bystep[label]` are `[proposed, accepted]`.
struct _MoveStats
    logscales::Vector{Float64}
    counts::Dict{Symbol, Vector{Int}}
    bystep::Dict{Any, Vector{Int}}
end

_MoveStats(d, τ) =
    _MoveStats(fill(log(τ), _nscales(d)), Dict(:warmup => [0, 0], :sampling => [0, 0]), Dict{Any, Vector{Int}}())

_nscales(::NormalStep) = 1
_nscales(d::ComponentwiseNormalStep) = d.n
_nscales(d) = 0

# Warmup calls since the last logged table: accepted proposals per move, the calls, their
# wall time, and the sampler step of the first of them.
mutable struct _MoveTally
    accepted::Vector{Int}
    calls::Int
    seconds::Float64
    first::Union{Nothing, Int}
end
_MoveTally(n) = _MoveTally(zeros(Int, n), 0, 0.0, nothing)

"""
    MoveSet(post::VLBIPosterior, moves; space = nothing, θ0 = nothing, rounds = 1,
            initial_scale = 0.1, target_accept = 0.45, output = nothing,
            rng = Random.default_rng())

A set of [`AbstractMove`](@ref)s in the latent space `space` of `post`, with their step scales
and acceptance statistics.

  - `rounds`: proposals of each move per call.
  - `initial_scale`: the starting step scale of moves with a [`NormalStep`](@ref) or
    [`ComponentwiseNormalStep`](@ref); ignored for other moves, whose entries may be
    `nothing`.
  - `target_accept`: the acceptance rate the step scales adapt toward during warmup.
  - `θ0`: a parameter point. If given, construction checks there that every
    [`likelihood_invariant`](@ref) move leaves the likelihood unchanged, and that the
    components of every componentwise move are independent, with steps drawn from `rng`.
  - `output`: a file path. If given, [`move_summary`](@ref) is serialized there after every
    call.

`rounds`, `initial_scale` and `target_accept` take one number for all moves or a vector with
one value per move. A componentwise move must be [`likelihood_invariant`](@ref).

    ms(x, tpost, phase, rng; step = nothing, total = nothing, report = true) -> x′

Run the moves on the latent point `x` of `tpost`, which must use the same space, in the
sampler phase `phase` (`:warmup` or `:sampling`). Call it between the steps of any MCMC
sampler. Returns `x` itself if no proposal was accepted. A proposal with a NaN acceptance
ratio is rejected. During warmup the step scales adapt, and with `report` the call logs the
acceptance rates and step scales since the last log, labeled with `step` out of `total`.
Throws if an accepted
[`likelihood_invariant`](@ref) move changed the likelihood.

    ms(state, tpost, info::NamedTuple, rng) -> state

The `between_chunks` hook of [`ReactantNUTS`](@ref): moves `state.position`, reading `phase`,
`step`, `total` and `report` from `info`.

The manual page "Writing and using moves" describes how proposals are accepted and run on
the device.
"""
struct MoveSet{M <: Tuple, D <: Tuple, C}
    moves::M
    dists::D
    ctx::C
    rounds::Vector{Int}
    initial_scale::Vector{Float64}
    target_accept::Vector{Float64}
    stats::Vector{_MoveStats}
    free::Vector{Bool}
    output::Union{Nothing, String}
    compiled::Base.RefValue{Any}
    seconds::Base.RefValue{Float64}
    tally::_MoveTally
end

_per_move(v::Real, n, name) = fill(v, n)
function _per_move(v::AbstractVector, n, name)
    length(v) == n || throw(ArgumentError("$name has $(length(v)) entries for $n moves"))
    return collect(v)
end

Base.@nospecializeinfer function MoveSet(
        @nospecialize(post::VLBIPosterior), @nospecialize(moves); space = nothing, θ0 = nothing, rounds = 1,
        initial_scale = 0.1, target_accept = 0.45, output = nothing,
        rng::AbstractRNG = Random.default_rng()
    )
    moves = Tuple(moves)
    isempty(moves) && throw(ArgumentError("MoveSet needs at least one move"))
    all(m -> m isa AbstractMove, moves) ||
        throw(ArgumentError("every move must be an AbstractMove"))
    names = map(move_name, moves)
    allunique(names) || throw(ArgumentError("move names must be distinct, got $(collect(names))"))
    dists = map(_step_distribution, moves)
    for (m, d) in zip(moves, dists)
        (traceable(m) && !(d isa NormalStep)) && throw(
            ArgumentError("move $(move_name(m)) is traceable but does not have a NormalStep")
        )
    end
    n = length(moves)
    r = _per_move(rounds, n, "rounds")
    all(>=(0), r) && any(>(0), r) ||
        throw(ArgumentError("rounds must be non-negative with at least one positive, got $r"))
    ta = _per_move(target_accept, n, "target_accept")
    all(a -> 0 < a < 1, ta) ||
        throw(ArgumentError("target_accept must lie in (0, 1), got $ta"))
    τ0 = _per_move(initial_scale, n, "initial_scale")
    scales = map(1:n) do j
        _nscales(dists[j]) > 0 || return NaN
        (τ0[j] isa Real && τ0[j] > 0) || throw(
            ArgumentError("the initial_scale of move $(names[j]) must be positive, got $(τ0[j])")
        )
        return Float64(τ0[j])
    end
    view = CoordinateView(post, space)
    ctx = move_context(view, moves)
    free = [likelihood_invariant(m) for m in moves]
    for j in 1:n
        (dists[j] isa ComponentwiseNormalStep && !free[j]) && throw(
            ArgumentError(
                "move $(names[j]) is componentwise, which needs a move that leaves the " *
                    "likelihood unchanged"
            )
        )
    end
    (any(free) && space isa PT.StdNormal) && _check_std_density(view)
    if !isnothing(θ0)
        x0 = inverse(view.tbase, θ0)
        for j in 1:n
            step = _draw_step(dists[j], rng, scales[j])
            likelihood_invariant(moves[j]) && _check_invariance(moves[j], ctx, x0; step)
            dists[j] isa ComponentwiseNormalStep && _check_components(moves[j], x0, step, ctx)
        end
    end
    return MoveSet(
        moves, dists, ctx, Int.(r), scales, Float64.(ta),
        [_MoveStats(dists[j], scales[j]) for j in 1:n], free,
        isnothing(output) ? nothing : String(output), Ref{Any}(nothing), Ref(0.0), _MoveTally(n)
    )
end

# The change of the latent prior log density of `tb` (in the latent space `space`) from `x`
# to `x′`: `log N(x; 0, I)` in the StdNormal space, which `_check_std_density` verifies.
_prior_delta(::PT.StdNormal, tb, x, x′) = -sum((x′ .- x) .* (x′ .+ x)) / 2
_prior_delta(space, tb, x, x′) = _latent_logprior(tb, x′) - _latent_logprior(tb, x)
_latent_logprior(tb, x) = last(latent_pfwd_and_logdensity(tb.transform, x))

# Throws unless the StdNormal log density of `view` is the log-likelihood plus
# `log N(x; 0, I)`, up to a constant.
function _check_std_density(view::CoordinateView; rtol::Real = 1.0e-8)
    tb = view.tbase
    rng = Random.Xoshiro(2)
    d(x) = logdensityof(tb, x) - loglikelihood(tb.lpost, transform(tb, x)) + sum(abs2, x) / 2
    d1, d2 = d(randn(rng, dimension(tb)) ./ 2), d(randn(rng, dimension(tb)) ./ 2)
    abs(d1 - d2) <= rtol * (abs(d1) + 1) || error(
        "the StdNormal log density is not the log-likelihood plus log N(x; 0, I): their " *
            "difference changes by $(abs(d1 - d2)) between two points, so invariant moves " *
            "cannot be accepted without the likelihood"
    )
    return nothing
end

# The change in the log-likelihood when `move` proposes from `x` with `step` (by default drawn
# from `rng` at scale 0.1); throws if it exceeds `rtol` relative to the log-likelihood.
Base.@nospecializeinfer function _check_invariance(
        @nospecialize(move::AbstractMove), @nospecialize(ctx), x;
        rng::AbstractRNG = Random.default_rng(), step = _draw_step(_step_distribution(move), rng, 0.1),
        rtol::Real = 1.0e-8
    )
    tb = ctx.view.tbase
    x′, _ = propose(move, x, step, ctx)
    l0 = loglikelihood(tb.lpost, transform(tb, x))
    l1 = loglikelihood(tb.lpost, transform(tb, x′))
    abs(l1 - l0) <= rtol * (abs(l0) + 1) || error(
        "move $(move_name(move)) changed the log-likelihood from $l0 to $l1; the model is " *
            "not invariant under it"
    )
    return abs(l1 - l0)
end

# The move context with its view over `tpost` and the arrays `data`.
_bind_context(data, tpost) = merge((; view = CoordinateView(tpost, PT.transport_node(tpost.transform))), data)

# What a `MoveSet` keeps between calls for positions like `position`: the moves' context
# arrays, and for a device position (Reactant extension) their device copies and the compiled
# programs.
_compile_moves(ms::MoveSet, tpost, position) = (; data = Base.structdiff(ms.ctx, NamedTuple{(:view,)}))

# The operations a `MoveSet` call runs through, around `tpost` and the cached `c`, for a host
# position; the Reactant extension defines the same set for a device position.
# `fused(x, ℓ, order, steps, logu)` runs the traceable moves' steps of one call and returns
# `(x, ℓ, logα, accepted)`.
Base.@nospecializeinfer function _move_kernels(@nospecialize(ms::MoveSet), @nospecialize(c), @nospecialize(tpost), @nospecialize(position))
    ctx = _bind_context(c.data, tpost)
    host_propose = (j, x, step) -> propose(ms.moves[j], x, step, ctx)
    k = (;
        load = p -> collect(Float64, vec(Array(p))),
        logdensity = Base.Fix1(logdensityof, tpost),
        prior_delta = (x, x′) -> _prior_delta(space(ctx.view), ctx.view.tbase, x, x′),
        propose = host_propose,
        host_propose,
        context = ctx,
        download = copy,
        upload = (xh, x) -> xh,
        store = _position_like,
    )
    tj = _traced_moves(ms)
    return merge(k, (; fused = (args...) -> _fused_sequential(k, tj, ms.free, args...)))
end
_position_like(x, position) = copyto!(similar(position), x)

_traced_moves(ms::MoveSet) = [j for j in eachindex(ms.moves) if traceable(ms.moves[j])]

# One step of move `j` with `step` from `x`, whose log density is `ℓ`, and the uniform
# `logu`: the point and its log density after the step, `log α`, and whether it was accepted.
Base.@nospecializeinfer function _transition(@nospecialize(k), free, j, x, ℓ, step, logu)
    x′, logdet = k.propose(j, x, step)
    ℓ′ = free ? ℓ + k.prior_delta(x, x′) : k.logdensity(x′)
    logα = ℓ′ - ℓ + logdet
    accepted = !isnan(logα) && logu < logα
    return accepted ? x′ : x, accepted ? ℓ′ : ℓ, logα, accepted
end

# The traceable moves' steps of one call: step `s` is the step `steps[s]` of move
# `tj[order[s]]`, with the uniform `logu[s]`.
Base.@nospecializeinfer function _fused_sequential(@nospecialize(k), tj, free, x, ℓ, order, steps, logu)
    logα = similar(steps)
    accepted = falses(length(steps))
    for (s, c) in pairs(order)
        j = tj[c]
        x, ℓ, logα[s], accepted[s] = _transition(k, free[j], j, x, ℓ, steps[s], logu[s])
    end
    return x, ℓ, logα, accepted
end

# The cache is keyed by `objectid(tpost)` rather than holding `tpost`, so a `MoveSet` kept
# after sampling does not keep the sampler's posterior alive.
Base.@nospecializeinfer function _kernels_for(@nospecialize(ms::MoveSet), @nospecialize(tpost), @nospecialize(position))
    c = ms.compiled[]
    if isnothing(c) || c.key != objectid(tpost) || c.T !== typeof(position)
        c = (; key = objectid(tpost), T = typeof(position), compiled = _compile_moves(ms, tpost, position))
        ms.compiled[] = c
    end
    return _move_kernels(ms, c.compiled, tpost, position)
end

"""
    update_context!(ms::MoveSet, key::Symbol, value) -> MoveSet

Copy `value` into the move-context array `key` of `ms` (same size). The moves read the new
value from the next call on; a compiled device program is updated without recompiling.
"""
function update_context!(ms::MoveSet, key::Symbol, value)
    haskey(ms.ctx, key) || throw(ArgumentError("the move context has no entry $key"))
    dest = getproperty(ms.ctx, key)
    dest isa AbstractArray || throw(ArgumentError("the move-context entry $key is not an array"))
    size(dest) == size(value) ||
        throw(DimensionMismatch("$key has size $(size(dest)), the new value $(size(value))"))
    copyto!(dest, value)
    c = ms.compiled[]
    if !isnothing(c)
        dev = getproperty(c.compiled.data, key)
        dev === dest || copyto!(dev, dest)
    end
    return ms
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
        st.logscales .= _adapt.(st.logscales, α, target, sweeps)
    end
    return st
end

# One round of the componentwise move `j` on the host: a step for every component, each
# accepted in turn from the prior alone. Returns the point, its log density and the number of
# components accepted.
Base.@nospecializeinfer function _componentwise_round(@nospecialize(ms::MoveSet), j, @nospecialize(k), x, ℓ, phase, rng)
    m = ms.moves[j]
    st = ms.stats[j]
    xh = k.download(x)
    x′, logdet = k.host_propose(j, xh, _draw_step(ms.dists[j], rng, _step_scale(ms.dists[j], st)))
    cs = component_coords(m)
    view = ms.ctx.view
    α = zeros(length(cs))
    accepted = falses(length(cs))
    for (i, c) in pairs(cs)
        xc = copy(xh)
        xc[c] = x′[c]
        Δ = _prior_delta(space(view), view.tbase, xh, xc)
        a = Δ + logdet[i]
        accepted[i] = !isnan(a) && log(rand(rng)) < a
        α[i] = _accept_prob(a)
        if accepted[i]
            xh = xc
            ℓ += Δ
        end
    end
    _record_components!(st, phase, α, accepted, ms.target_accept[j], [component_label(m, i) for i in eachindex(cs)])
    return any(accepted) ? k.upload(xh, x) : x, ℓ, count(accepted)
end

function _record!(st::_MoveStats, phase::Symbol, logα, accepted::Bool, target, label)
    _count!(st.counts[phase], accepted)
    isnothing(label) || _count!(get!(() -> [0, 0], st.bystep, label), accepted)
    phase === :warmup &&
        (st.logscales .= _adapt.(st.logscales, _accept_prob(logα), target, st.counts[:warmup][1]))
    return st
end

_accept_prob(logα) = isnan(logα) ? 0.0 : min(1.0, exp(logα))
_adapt(logτ, α, target, n) = clamp(logτ + (α - target) / n^0.6, log(1.0e-8), log(10.0))

function (ms::MoveSet)(
        position, tpost, phase::Symbol, rng::AbstractRNG;
        step = nothing, total = nothing, report::Bool = true
    )
    return _call(ms, position, tpost, phase, rng, step, total, report)
end

# Compiled once for all moves and posteriors, whose types are large; the calls inside dispatch
# to code compiled for each type.
Base.@nospecializeinfer function _call(
        @nospecialize(ms::MoveSet), @nospecialize(position), @nospecialize(tpost), phase::Symbol,
        rng::AbstractRNG, step, total, report::Bool
    )
    want = space(ms.ctx.view)
    have = _node_space(PT.transport_node(tpost.transform))
    typeof(have) == typeof(want) || throw(
        ArgumentError(
            "the moves act in the $(_space_name(want)) space but the sampler samples the " *
                "$(_space_name(have)) space"
        )
    )
    phase in (:warmup, :sampling) || throw(ArgumentError("unknown sampler phase $(repr(phase))"))
    t0 = Base.time()
    k = _kernels_for(ms, tpost, position)
    x = k.load(position)
    ℓ = k.logdensity(x)
    isfinite(ℓ) || error("the log density at the sampler's position is $ℓ before the moves")
    nacc = zeros(Int, length(ms.moves))
    tj = _traced_moves(ms)
    if !isempty(tj)
        order, steps, logu = Int[], Float64[], Float64[]
        for r in 1:maximum(ms.rounds[tj]), (c, j) in pairs(tj)
            r <= ms.rounds[j] || continue
            push!(order, c)
            push!(steps, _draw_step(ms.dists[j], rng, _step_scale(ms.dists[j], ms.stats[j])))
            push!(logu, log(rand(rng)))
        end
        x, ℓ, logα, accepted = k.fused(x, ℓ, order, steps, logu)
        for (s, c) in pairs(order)
            j = tj[c]
            label = step_label(ms.moves[j], steps[s])
            _record!(ms.stats[j], phase, logα[s], accepted[s], ms.target_accept[j], label)
            nacc[j] += accepted[s]
        end
    end
    for r in 1:maximum(ms.rounds), (j, m) in enumerate(ms.moves)
        (r <= ms.rounds[j] && !(j in tj)) || continue
        d = ms.dists[j]
        st = ms.stats[j]
        if d isa ComponentwiseNormalStep
            x, ℓ, a = _componentwise_round(ms, j, k, x, ℓ, phase, rng)
            nacc[j] += a
            continue
        end
        s = _draw_step(d, rng, _step_scale(d, st))
        x, ℓ, logα, accepted = _transition(k, ms.free[j], j, x, ℓ, s, log(rand(rng)))
        _record!(st, phase, logα, accepted, ms.target_accept[j], step_label(m, s))
        nacc[j] += accepted
    end
    changed = [move_name(m) for (j, m) in enumerate(ms.moves) if ms.free[j] && nacc[j] > 0]
    if !isempty(changed)
        ℓfull = k.logdensity(x)
        # Tight: a likelihood change under this tolerance is accepted as if it were none.
        abs(ℓfull - ℓ) <= 1.0e-9 * (abs(ℓfull) + 1) || error(
            "the log density after the moves is $ℓfull, but accepting them without the " *
                "likelihood gave $ℓ: one of the accepted moves " *
                "($(join(changed, ", "))) " *
                "changed the likelihood"
        )
    end
    isnothing(ms.output) || serialize(ms.output, move_summary(ms))
    moved = any(>(0), nacc) ? k.store(x, position) : position
    dt = Base.time() - t0
    ms.seconds[] += dt
    phase === :warmup && _log_warmup_moves!(ms, nacc, dt, step, total, report)
    return moved
end

function (ms::MoveSet)(state, tpost, info::NamedTuple, rng::AbstractRNG)
    x = state.position
    x′ = ms(
        x, tpost, info.phase, rng;
        step = get(info, :step, nothing), total = get(info, :total, nothing),
        report = get(info, :report, true)
    )
    x′ === x || (state.position = x′)
    return state
end

function _log_warmup_moves!(ms::MoveSet, nacc, dt, step, total, report)
    t = ms.tally
    t.calls == 0 && (t.first = step)
    t.accepted .+= nacc
    t.calls += 1
    t.seconds += dt
    report || return
    at = isnothing(step) ? "warmup" : "warmup step $step" * (isnothing(total) ? "" : "/$total")
    calls = t.calls == 1 ? "1 call" :
        isnothing(t.first) ? "$(t.calls) calls" : "$(t.calls) calls from step $(t.first)"
    @info "moves after $at ($calls, $(@sprintf("%.2f", t.seconds)) s)\n" *
        _moves_table(ms, t.accepted, t.calls)
    fill!(t.accepted, 0)
    t.calls = 0
    t.seconds = 0.0
    return
end

# One row per move: proposals accepted over `calls` calls, proposals made, the acceptance
# rate, and the step scale (the geometric mean over components for a componentwise move;
# blank for moves without one).
function _moves_table(ms::MoveSet, nacc, calls)
    rows = map(enumerate(ms.moves)) do (j, m)
        d = ms.dists[j]
        st = ms.stats[j]
        nprop = calls * (d isa ComponentwiseNormalStep ? ms.rounds[j] * d.n : ms.rounds[j])
        τ = isempty(st.logscales) ? nothing : exp(sum(st.logscales) / length(st.logscales))
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

Statistics for each move: `name`; the current step scale `τ` (a vector for a componentwise
move, `nothing` for a move without one); `rounds`; the `proposed` and `accepted` counts in
`warmup` and `sampling`; and `bystep`, the counts per [`step_label`](@ref) (per
[`component_label`](@ref) for a componentwise move) as `label => [proposed, accepted]`.
"""
function move_summary(ms::MoveSet)
    return map(enumerate(ms.moves)) do (j, m)
        st = ms.stats[j]
        (;
            name = move_name(m),
            τ = _step_scale(ms.dists[j], st),
            rounds = ms.rounds[j],
            warmup = (proposed = st.counts[:warmup][1], accepted = st.counts[:warmup][2]),
            sampling = (proposed = st.counts[:sampling][1], accepted = st.counts[:sampling][2]),
            bystep = copy(st.bystep),
        )
    end
end

_step_scale(::NormalStep, st) = exp(only(st.logscales))
_step_scale(::ComponentwiseNormalStep, st) = exp.(st.logscales)
_step_scale(d, st) = nothing

# `post` with data simulated at the parameters `θ`.
function _with_data(rng, post, θ)
    data = simulate_observation(rng, post, θ)
    return setproperties(post, (; data, lklhds = map(makelikelihood, data)))
end

# The KS statistics and their critical value of Geweke's joint-distribution test of the
# kernel of one call of `ms` (see `check_move`) on the latent coordinates `idx`, the
# projections onto the vectors `dirs` and the log-likelihood, and the acceptance rate of
# each move; an error if a statistic exceeds the critical value.
Base.@nospecializeinfer function _check_geweke(@nospecialize(ms::MoveSet), idx, dirs, n, steps, ks_alpha, rng)
    tbase = ms.ctx.view.tbase
    post = tbase.lpost
    isnothing(post.dataimg) || throw(
        ArgumentError("the joint-distribution test cannot simulate the image data of the posterior")
    )
    moved = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        θ = prior_sample(rng, post)
        x = inverse(tbase, θ)
        py = _with_data(rng, post, θ)
        for t in 1:steps
            t > 1 && (py = _with_data(rng, post, transform(tbase, x)))
            x = ms(x, TransformedVLBIPosterior(py, tbase.transform), :sampling, rng)
        end
        moved[i] = vcat(x[idx], dot.(dirs, Ref(x)), loglikelihood(py, transform(tbase, x)))
    end
    fresh = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        θ = prior_sample(rng, post)
        x = inverse(tbase, θ)
        fresh[i] = vcat(x[idx], dot.(dirs, Ref(x)), loglikelihood(_with_data(rng, post, θ), θ))
    end
    ks = [_ks_statistic(getindex.(moved, k), getindex.(fresh, k)) for k in eachindex(first(fresh))]
    crit = _ks_critical(ks_alpha, length(ks), n)
    all(<=(crit), ks) || error(
        "the kernel of $(join(map(move_name, ms.moves), ", ")) does not leave the posterior " *
            "invariant: KS statistics $(round.(ks; digits = 3)) on coordinates $idx, " *
            "$(length(dirs)) move directions and the log-likelihood against a critical " *
            "value $(round(crit; digits = 3))"
    )
    return ks, crit, [st.counts[:sampling][2] / st.counts[:sampling][1] for st in ms.stats]
end

"""
    check_moveset(ms::MoveSet; ngeweke, geweke_steps = 10, rng = Random.default_rng(),
                  ncoords = 5, ks_alpha = 1e-3) -> NamedTuple

Geweke's joint-distribution test (see [`check_move`](@ref)) of one call of `ms` at its current
step scales, with `ngeweke` replicates of `geweke_steps` steps. It compares up to `ncoords`
of the coordinates each move changes, and the position along each move's direction. Returns
`geweke`, `geweke_critical` and `geweke_accept` (one rate per move). Throws if a KS statistic
exceeds its critical value.
"""
Base.@nospecializeinfer function check_moveset(
        @nospecialize(ms::MoveSet); ngeweke::Integer, geweke_steps::Integer = 10,
        rng::AbstractRNG = Random.default_rng(), ncoords::Integer = 5, ks_alpha::Real = 1.0e-3
    )
    ngeweke > 0 || throw(ArgumentError("ngeweke must be positive, got $ngeweke"))
    view = ms.ctx.view
    post = view.tbase.lpost
    run = MoveSet(
        post, ms.moves; space = space(view), rounds = ms.rounds,
        initial_scale = ms.initial_scale, target_accept = ms.target_accept
    )
    for (a, b) in zip(run.stats, ms.stats)
        a.logscales .= b.logscales
    end
    xs = prior_sample(rng, view.tbase, 10)
    idxs, dirs = Vector{Int}[], Vector{Float64}[]
    for j in eachindex(ms.moves)
        d = ms.dists[j]
        changed, direction = _first_change(ms.moves[j], d, xs, ms.ctx, rng, _step_scale(d, ms.stats[j]))
        push!(idxs, _spread(changed, ncoords))
        push!(dirs, direction)
    end
    idx = sort!(unique(reduce(vcat, idxs)))
    geweke, geweke_critical, geweke_accept =
        _check_geweke(run, idx, dirs, ngeweke, geweke_steps, ks_alpha, rng)
    return (; geweke, geweke_critical, geweke_accept)
end
