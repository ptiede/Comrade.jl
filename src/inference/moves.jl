# Metropolis–Hastings moves in the base latent coordinates of a posterior (flat or
# StdNormal), the parameter view they read and write through, and checks of their
# correctness.

export CoordinateView, AbstractMove, NormalStep, ComponentwiseNormalStep, SymmetricDiscreteStep,
    CompensatedMove, check_move

"""
    CoordinateView(post::VLBIPosterior, space = nothing)

The parameters of `post` in a latent space: the flat space of [`asflat`](@ref)
(`space = nothing`) or `ProbabilityTransports.StdNormal()`. A parameter is named by its path,
such as `(:sky, :σa)` or `(:instrument, :lg1)`, and is read and written with
[`coords`](@ref), [`value`](@ref) and [`latent`](@ref). `view.tbase` is the transformed
posterior.
"""
struct CoordinateView{T, N}
    tbase::T
    root::N
end

function CoordinateView(post::VLBIPosterior, space = nothing)
    space isa Union{Nothing, PT.StdNormal} || throw(
        ArgumentError("moves act in the flat or StdNormal latent space, not $(typeof(space))")
    )
    tbase = isnothing(space) ? asflat(post) : transport_to(post, space)
    return CoordinateView(tbase, PT.transport_node(tbase.transform))
end

space(view::CoordinateView) = _node_space(view.root)
_node_space(root) = root isa PT.AbstractTransport ? PT.space(root) : nothing
_space_name(::Nothing) = "flat"
_space_name(s) = string(nameof(typeof(s)))

Base.show(io::IO, v::CoordinateView) =
    print(io, "CoordinateView($(_space_name(space(v))), dimension $(dimension(v.tbase)))")

_children(t::TV.TransformTuple) = getfield(t, :inner)
_children(t::PT.TupleTransport) = getfield(t, :transports)
_children(t) = nothing
_node_dimension(t::PT.AbstractTransport) = PT.dimension(t)
_node_dimension(t) = TV.dimension(t)

# The node at `path` below `t` and the latent coordinates it consumes, for a `t` whose
# coordinates start after `offset`.
function _locate_node(t, path::Tuple, offset::Int, full::Tuple)
    isempty(path) && return t, (offset + 1):(offset + _node_dimension(t))
    children = _children(t)
    isnothing(children) && throw(
        ArgumentError("no parameter at $(full): $(full[1:(length(full) - length(path))]) is a $(nameof(typeof(t))) leaf")
    )
    for k in keys(children)
        k == first(path) && return _locate_node(children[k], Base.tail(path), offset, full)
        offset += _node_dimension(children[k])
    end
    throw(
        ArgumentError(
            "no parameter $(first(path)) at $(full[1:(length(full) - length(path))]); " *
                "found $(collect(keys(children)))"
        )
    )
end

"""
    node(view::CoordinateView, path) -> transform node

The transform node of the parameter at `path`.
"""
node(view::CoordinateView, path::Tuple) = first(_locate_node(view.root, path, 0, path))

"""
    coords(view::CoordinateView, path) -> UnitRange

The latent coordinates of the parameter at `path`.
"""
coords(view::CoordinateView, path::Tuple) = last(_locate_node(view.root, path, 0, path))

_block_value(t::PT.AbstractTransport, y) = PT.latent_pfwd(t, y)
_block_value(t::TV.ScalarTransform, y) = TV.transform(t, _rget(y, firstindex(y)))
_block_value(t, y) = TV.transform(t, y)
_block_latent(t::PT.AbstractTransport, v) = PT.latent_pback(t, v)
_block_latent(t, v) = TV.inverse(t, v)

"""
    value(view::CoordinateView, x, path)

The value of the parameter at `path` at the latent point `x`.
"""
function value(view::CoordinateView, x, path::Tuple)
    t, r = _locate_node(view.root, path, 0, path)
    return _block_value(t, x[r])
end

"""
    latent(view::CoordinateView, path, v) -> Vector

The latent coordinates of the parameter at `path` that give it the value `v`.
"""
latent(view::CoordinateView, path::Tuple, v) = _latent_vector(_block_latent(node(view, path), v))

_latent_vector(y::Number) = [y]
_latent_vector(y) = vec(y)

# `latent` of the block node `t`, allocated like `x` so a traced `x` gives a traced block.
_block_latent_like(t::PT.AbstractTransport, v, x) = PT.latent_pback!(similar(x, PT.dimension(t)), t, v)
_block_latent_like(t, v, x) = _block_latent(t, v)
_block_latent_like(::TV.ArrayTransformation{TV.Identity}, v, x) = vec(v)

# --- the move protocol ------------------------------------------------------------------

"""
    AbstractMove

Supertype of Metropolis–Hastings moves on the latent point of a posterior. A move implements
[`propose`](@ref), [`step_distribution`](@ref) and [`move_name`](@ref), and may implement
[`likelihood_invariant`](@ref), [`step_label`](@ref), [`component_label`](@ref),
[`context_data`](@ref) and [`traceable`](@ref).

`propose` must not modify `x`. Proposing with the reverse step from the proposal must return
to `x` with the opposite log-determinant.
[`check_move`](@ref) tests a move.
"""
abstract type AbstractMove end

"""
    NormalStep()

A scalar step `u ~ N(0, τ²)`, where `τ` is the move's step scale. The reverse step is `−u`.
"""
struct NormalStep end

"""
    ComponentwiseNormalStep(n)

A vector step of `n` components, `uᵢ ~ N(0, τᵢ²)`, each with its own step scale. The reverse
step is `−u`. Each component is accepted or rejected on its own, so the move must implement
[`component_coords`](@ref), the components must change disjoint coordinates, and component
`i` of the proposal must depend only on `uᵢ` and the coordinates of component `i`.
[`propose`](@ref) returns one log-determinant per component.
"""
struct ComponentwiseNormalStep
    n::Int
    function ComponentwiseNormalStep(n)
        n > 0 || throw(ArgumentError("a componentwise move needs at least one component, got $n"))
        return new(n)
    end
end

"""
    SymmetricDiscreteStep(draw, reverse)

A step drawn by `draw(rng)`. `reverse(step)` is the reverse step, which `draw` must return
with the same probability as `step`.
"""
struct SymmetricDiscreteStep{D, R}
    draw::D
    reverse::R
end

"""
    propose(move, x, step, ctx)

The proposal of `move` from the latent point `x`: `(x′, logdet)` with
`logdet = log|det ∂x′/∂x|` (one per component for a [`ComponentwiseNormalStep`](@ref)).
`ctx` comes from [`move_context`](@ref).
"""
function propose end

"""
    step_distribution(move)

The step distribution of `move`: a [`NormalStep`](@ref), [`ComponentwiseNormalStep`](@ref)
or [`SymmetricDiscreteStep`](@ref).
"""
function step_distribution end

"""
    move_name(move) -> String

The name of `move` in logs and statistics.
"""
function move_name end

"""
    likelihood_invariant(move) -> Bool

Whether `move` leaves the likelihood unchanged (default `false`). [`MoveSet`](@ref) accepts
such a move from the change in the prior alone, without evaluating the likelihood, then
checks that the likelihood did not change.
"""
likelihood_invariant(::AbstractMove) = false

"""
    context_data(move, view) -> NamedTuple

Arrays that `move` reads from `ctx` in [`propose`](@ref) (default none).
[`update_context!`](@ref) replaces their contents.
"""
context_data(::AbstractMove, view) = NamedTuple()

"""
    traceable(move) -> Bool

Whether [`propose`](@ref) for `move` compiles with Reactant (default `false`). On a device
position, [`MoveSet`](@ref) runs such a proposal on the device, with `ctx.view` over the
device posterior. A traceable move has a [`NormalStep`](@ref) and reads every array it needs
from `ctx` rather than capturing it.
"""
traceable(::AbstractMove) = false

"""
    component_coords(move) -> Vector{Vector{Int}}

The latent coordinates each component of a [`ComponentwiseNormalStep`](@ref) move changes.
"""
function component_coords end

"""
    component_label(move, i)

The name of component `i` of a [`ComponentwiseNormalStep`](@ref) move in its statistics
(default `i`).
"""
component_label(::AbstractMove, i) = i

"""
    step_label(move, step)

A key under which [`move_summary`](@ref) counts the proposals and acceptances of `step`, or
`nothing` (default) to count only the move's totals.
"""
step_label(::AbstractMove, step) = nothing

function _step_distribution(m::AbstractMove)
    d = step_distribution(m)
    d isa Union{NormalStep, ComponentwiseNormalStep, SymmetricDiscreteStep} || throw(
        ArgumentError(
            "the step distribution of move $(move_name(m)) is a $(typeof(d)), not a " *
                "NormalStep, ComponentwiseNormalStep or SymmetricDiscreteStep"
        )
    )
    return d
end

_draw_step(::NormalStep, rng, τ) = τ * randn(rng)
_draw_step(d::ComponentwiseNormalStep, rng, τ) = τ .* randn(rng, d.n)
_draw_step(d::SymmetricDiscreteStep, rng, τ) = d.draw(rng)

_reverse_step(::Union{NormalStep, ComponentwiseNormalStep}, u) = -u
_reverse_step(d::SymmetricDiscreteStep, step) = d.reverse(step)

"""
    move_context(view::CoordinateView, moves) -> NamedTuple

The `ctx` argument of [`propose`](@ref): `(; view, ...)` merged with the `context_data` of
every move.
"""
function move_context(view::CoordinateView, moves)
    return foldl((c, m) -> merge(c, context_data(m, view)), moves; init = (; view))
end

_rget(x, i) = ComradeBase.rgetindex(x, i)

function _with_coordinate(x, i::Int, v)
    x′ = copy(x)
    ComradeBase.rsetindex!(x′, v, i)
    return x′
end

"""
    CompensatedMove(name, view, shift, block, compensate; index = 1,
                    logdet = (vb, x, x′, ctx) -> 0.0, likelihood_invariant = false,
                    context = (;), traceable = false)

A [`NormalStep`](@ref) move that adds the step to the `index`-th latent coordinate of the
parameter at path `shift`, and sets the parameter at path `block` to
`compensate(vb, x, x′, ctx)`. Here `vb` is the value of the block at `x`, and `x′` is `x` with
only the shift applied. `compensate` keeps fixed something the likelihood depends on, such as
a product of parameters.

  - `logdet(vb, x, x′, ctx)`: `log|det ∂x′/∂x|` of the whole move, given the complete
    proposal `x′`.
  - `likelihood_invariant`: see [`likelihood_invariant`](@ref).
  - `context`: a `NamedTuple` of arrays added to `ctx`.
  - `traceable`: see [`traceable`](@ref). `compensate` and `logdet` must then read
    `ctx.view` and the `context` arrays from `ctx`.

Throws if `shift` or `block` is not in `view`, or if the shifted coordinate lies in `block`.
"""
struct CompensatedMove{B <: Tuple, C, L, X <: NamedTuple} <: AbstractMove
    name::String
    ishift::Int
    block::B
    compensate::C
    logdet::L
    likelihood_invariant::Bool
    context::X
    traceable::Bool
end

function CompensatedMove(
        name::AbstractString, view::CoordinateView, shift::Tuple, block::Tuple, compensate;
        index::Integer = 1, logdet = (vb, x, x′, ctx) -> 0.0, likelihood_invariant::Bool = false,
        context::NamedTuple = (;), traceable::Bool = false
    )
    rs = coords(view, shift)
    1 <= index <= length(rs) || throw(
        ArgumentError("index $index is outside the $(length(rs)) latent coordinates of $shift")
    )
    rb = coords(view, block)
    i = rs[index]
    i in rb && throw(ArgumentError("the shifted coordinate of $shift lies in the block $block"))
    return CompensatedMove(String(name), i, block, compensate, logdet, likelihood_invariant, context, traceable)
end

step_distribution(::CompensatedMove) = NormalStep()
move_name(m::CompensatedMove) = m.name
likelihood_invariant(m::CompensatedMove) = m.likelihood_invariant
context_data(m::CompensatedMove, view) = m.context
traceable(m::CompensatedMove) = m.traceable

function propose(m::CompensatedMove, x, u, ctx)
    x′ = _with_coordinate(x, m.ishift, _rget(x, m.ishift) + u)
    t, rb = _locate_node(ctx.view.root, m.block, 0, m.block)
    vb = _block_value(t, x[rb])
    x′[rb] = _latent_vector(_block_latent_like(t, m.compensate(vb, x, x′, ctx), x′))
    return x′, m.logdet(vb, x, x′, ctx)
end

# --- checks -----------------------------------------------------------------------------

# Two-sample Kolmogorov–Smirnov statistic.
function _ks_statistic(a::AbstractVector, b::AbstractVector)
    sa, sb = sort(a), sort(b)
    d = 0.0
    for v in vcat(sa, sb)
        fa = searchsortedlast(sa, v) / length(sa)
        fb = searchsortedlast(sb, v) / length(sb)
        d = max(d, abs(fa - fb))
    end
    return d
end

# Throws unless the components of a componentwise move are independent at `x`: disjoint
# coordinates, nothing changed outside them, and component `i` of `step` alone changing only
# component `i`'s coordinates, as in the full step and with the same log-determinant.
function _check_components(move::AbstractMove, x, step, ctx; rtol::Real = 1.0e-8)
    name = move_name(move)
    cs = component_coords(move)
    length(cs) == length(step) ||
        error("move $name has $(length(cs)) components but a step of length $(length(step))")
    every = reduce(vcat, cs; init = Int[])
    allunique(every) || error("the components of move $name share latent coordinates")
    x′, ld = propose(move, x, step, ctx)
    tol = rtol * (maximum(abs, x; init = 0.0) + 1)
    outside = setdiff(eachindex(x), every)
    maximum(abs, x′[outside] .- x[outside]; init = 0.0) <= tol ||
        error("move $name changes latent coordinates outside its components")
    for i in eachindex(cs)
        u = zero(step)
        u[i] = step[i]
        xi, ldi = propose(move, x, u, ctx)
        rest = setdiff(every, cs[i])
        (
            maximum(abs, xi[cs[i]] .- x′[cs[i]]; init = 0.0) <= tol &&
                maximum(abs, xi[rest] .- x[rest]; init = 0.0) <= tol &&
                isapprox(ldi[i], ld[i]; atol = rtol * (abs(ld[i]) + 1))
        ) || error(
            "component $(component_label(move, i)) of move $name depends on the other " *
                "components' steps or changes their coordinates"
        )
    end
    return nothing
end

# Up to `n` entries of `c`, evenly spaced.
_spread(c, n) = c[unique(round.(Int, range(1, length(c), length = min(n, length(c)))))]

# The critical value of `m` two-sample KS statistics of samples of size `n` at total level
# `alpha` (Bonferroni).
_ks_critical(alpha, m, n) = sqrt(-log(alpha / m / 2) / 2) * sqrt(2 / n)

_logprior_latent(tbase, post, x) =
    logdensityof(tbase, x) - loglikelihood(post, transform(tbase, x))

# The coordinates changed by the first of up to 100 steps of `move` (from the points of `xs`
# in turn) that changes its point, and the unit direction of the change.
function _first_change(move, d, xs, ctx, rng, τ)
    for i in 1:100
        x = xs[mod1(i, length(xs))]
        step = _draw_step(d, rng, τ)
        x′, _ = propose(move, x, step, ctx)
        changed = findall(x′ .!= x)
        isempty(changed) || return changed, normalize(x′ .- x)
    end
    error("100 steps of move $(move_name(move)) left the point unchanged")
end

"""
    check_move(move, post, θs; space = nothing, rng = Random.default_rng(), τ = 0.1,
               rtol = 1e-8, h = 1e-6, logdet_atol = 1e-5, nprior = 0, nsteps = 5,
               ncoords = 5, ks_alpha = 1e-3, ngeweke = 0, geweke_steps = 10) -> NamedTuple

Test `move` on `post` in the latent space `space`, at each parameter point in `θs`, with steps
at scale `τ`. Throws an error naming the first failed check:

  - reversal: the reverse step from the proposal returns to the start, with the opposite
    log-determinant;
  - log-determinant: `logdet` matches a finite-difference Jacobian (step `h`);
  - invariance: if [`likelihood_invariant`](@ref)`(move)`, the log-likelihood is unchanged;
  - components: for a [`ComponentwiseNormalStep`](@ref) move, the components are independent:
    disjoint coordinates, and component `i` of the step changes only component `i`.

Two optional tests compare samples with two-sample Kolmogorov–Smirnov tests at level
`ks_alpha`, on up to `ncoords` of the coordinates the move changes:

  - `nprior > 0`: `nprior` prior draws each take `nsteps` steps of the move targeting the
    prior, and are compared with fresh prior draws.
  - `ngeweke > 0`: Geweke's (2004) joint-distribution test of the kernel [`MoveSet`](@ref)
    runs. Each of `ngeweke` replicates draws parameters from the prior and data from the model
    ([`simulate_observation`](@ref)), then alternates `geweke_steps` move steps with new data.
    The results are compared with fresh draws of parameters and data, also along the move's
    direction and on the log-likelihood. Arrays in the move context are not rebuilt for the
    new data.

Returns the largest reversal error, log-determinant error and log-likelihood change; the KS
statistics `ks` and `geweke` with their critical values; and `geweke_accept`, the acceptance
rate in the Geweke test. A move that rarely moves passes the Geweke test untested, so check
that `geweke_accept` is not near zero.
"""
function check_move(
        move::AbstractMove, post::VLBIPosterior, θs; space = nothing,
        rng::AbstractRNG = Random.default_rng(), τ::Real = 0.1, rtol::Real = 1.0e-8,
        h::Real = 1.0e-6, logdet_atol::Real = 1.0e-5, nprior::Integer = 0, nsteps::Integer = 5,
        ncoords::Integer = 5, ks_alpha::Real = 1.0e-3, ngeweke::Integer = 0,
        geweke_steps::Integer = 10
    )
    view = CoordinateView(post, space)
    ctx = move_context(view, (move,))
    d = _step_distribution(move)
    name = move_name(move)
    tbase = view.tbase
    rev, lderr, dll = 0.0, 0.0, 0.0
    xs = [inverse(tbase, θ) for θ in θs]
    for x in xs
        step = _draw_step(d, rng, τ)
        x′, ld = propose(move, x, step, ctx)
        x″, ld′ = propose(move, x′, _reverse_step(d, step), ctx)
        e = maximum(abs, x″ .- x) / (maximum(abs, x) + 1)
        rev = max(rev, e)
        d isa ComponentwiseNormalStep && _check_components(move, x, step, ctx; rtol)
        ld, ld′ = sum(ld), sum(ld′)
        (e <= rtol && isapprox(ld′, -ld; atol = logdet_atol)) || error(
            "move $name is not reversed by its reverse step at step $step: |Δx| = $e, " *
                "logdet $ld then $ld′ (should sum to 0)"
        )
        T = findall(x′ .!= x)
        isempty(T) && error("move $name left the point unchanged at step $step")
        J = similar(x, length(T), length(T))
        for (c, j) in enumerate(T)
            xp, xm = copy(x), copy(x)
            xp[j] += h
            xm[j] -= h
            J[:, c] = (propose(move, xp, step, ctx)[1][T] .- propose(move, xm, step, ctx)[1][T]) ./ (2h)
        end
        ldfd = logabsdet(J)[1]
        lderr = max(lderr, abs(ldfd - ld))
        abs(ldfd - ld) <= logdet_atol * max(1, abs(ld)) || error(
            "move $name reports logdet = $ld, but the finite-difference Jacobian over its " *
                "$(length(T)) changed coordinates gives $ldfd"
        )
        likelihood_invariant(move) &&
            (dll = max(dll, _check_invariance(move, ctx, x; step, rtol)))
    end
    if nprior > 0 || ngeweke > 0
        changed, direction = _first_change(move, d, xs, ctx, rng, τ)
    end
    ks, crit = Float64[], NaN
    if nprior > 0
        idx = _spread(changed, ncoords)
        crit = _ks_critical(ks_alpha, length(idx), nprior)
        draw() = inverse(tbase, prior_sample(rng, post))
        moved = map(1:nprior) do _
            x = draw()
            ℓ = _logprior_latent(tbase, post, x)
            for _ in 1:nsteps
                x′, ld = propose(move, x, _draw_step(d, rng, τ), ctx)
                ℓ′ = _logprior_latent(tbase, post, x′)
                if log(rand(rng)) < ℓ′ - ℓ + sum(ld)
                    x, ℓ = x′, ℓ′
                end
            end
            x[idx]
        end
        fresh = [draw()[idx] for _ in 1:nprior]
        ks = [_ks_statistic(getindex.(moved, k), getindex.(fresh, k)) for k in eachindex(idx)]
        all(<=(crit), ks) || error(
            "move $name does not leave the prior invariant: KS statistics $(round.(ks; digits = 3)) " *
                "on coordinates $idx against a critical value $(round(crit; digits = 3))"
        )
    end
    gks, gcrit, gacc = Float64[], NaN, NaN
    if ngeweke > 0
        ms = MoveSet(post, (move,); space, initial_scale = τ)
        gks, gcrit, gacc = _check_geweke(
            ms, _spread(changed, ncoords), [direction], ngeweke, geweke_steps, ks_alpha, rng
        )
        gacc = only(gacc)
    end
    return (;
        reversal = rev, logdet = lderr, loglikelihood = dll, ks, ks_critical = crit,
        geweke = gks, geweke_critical = gcrit, geweke_accept = gacc,
    )
end
