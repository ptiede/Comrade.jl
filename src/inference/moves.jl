# Metropolis–Hastings moves in the base latent coordinates of a posterior (flat or
# StdNormal), the parameter view they read and write through, and checks of their
# correctness.

export CoordinateView, AbstractMove, RandomWalk, DiscreteSymmetric, CompensatedMove, check_move

"""
    CoordinateView(post::VLBIPosterior, space = nothing)

The parameters of `post` addressed in the base latent space `space`: `nothing` for the
flat space ([`asflat`](@ref)) or `ProbabilityTransports.StdNormal()`. A parameter is named
by its path, a tuple of field names such as `(:sky, :σa)` or `(:instrument, :lg1)`.

  - [`coords`](@ref)`(view, path)`: the latent coordinates of the parameter;
  - [`value`](@ref)`(view, x, path)`: its value at the latent point `x`;
  - [`latent`](@ref)`(view, path, v)`: the latent coordinates of the value `v`.

`view.tbase` is the transformed posterior of the base space.
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

space(view::CoordinateView) = view.root isa PT.AbstractTransport ? PT.space(view.root) : nothing

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
_block_value(t::TV.ScalarTransform, y) = TV.transform(t, only(y))
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

# --- the move protocol ------------------------------------------------------------------

"""
    AbstractMove

A Metropolis–Hastings move of the base latent coordinates of a posterior. A move implements

  - [`propose`](@ref)`(move, x, step, ctx) -> (x′, logdet)`: the moved point and
    `log|det ∂x′/∂x|`, with `ctx` from [`move_context`](@ref);
  - `step_kind(move)`: [`RandomWalk`](@ref) (a scalar step `u ~ N(0, τ²)`) or
    [`DiscreteSymmetric`](@ref);
  - `move_name(move)`.

and optionally `draw_step(move, rng, τ)` and `reverse_step(move, step)` (required for
`DiscreteSymmetric`), `is_invariant(move)` (default `true`: the likelihood does not change)
and `context_data(move, view)` (arrays the move reads from `ctx`).

A proposal and the proposal with `reverse_step` must invert each other, and the step
distribution must be symmetric under `reverse_step`, so that
`log α = ℓ(x′) − ℓ(x) + logdet` with `ℓ` the log density in the base coordinates.
[`check_move`](@ref) tests these properties.
"""
abstract type AbstractMove end

"""
    RandomWalk(initial_scale)

The step kind of a move whose step is a scalar `u ~ N(0, τ²)`, with `τ` starting at
`initial_scale`; `u` and `−u` are reverse steps.
"""
struct RandomWalk
    initial_scale::Float64
    function RandomWalk(initial_scale)
        initial_scale > 0 ||
            throw(ArgumentError("initial_scale must be positive, got $initial_scale"))
        return new(initial_scale)
    end
end

"""
    DiscreteSymmetric()

The step kind of a move that draws its step from a finite set with `draw_step(move, rng, τ)`
(`τ` unused), where a step and `reverse_step(move, step)` are equally likely.
"""
struct DiscreteSymmetric end

"""
    propose(move, x, step, ctx) -> (x′, logdet)

The point `move` takes the base latent point `x` to with `step`, and `log|det ∂x′/∂x|`.
"""
function propose end

function step_kind end
function move_name end
is_invariant(::AbstractMove) = true
context_data(::AbstractMove, view) = NamedTuple()

draw_step(m::AbstractMove, rng::AbstractRNG, τ) = _draw_step(step_kind(m), m, rng, τ)
_draw_step(::RandomWalk, m, rng, τ) = τ * randn(rng)
_draw_step(::DiscreteSymmetric, m, rng, τ) =
    throw(ArgumentError("move $(move_name(m)) has discrete steps but no draw_step method"))

reverse_step(m::AbstractMove, step) = _reverse_step(step_kind(m), m, step)
_reverse_step(::RandomWalk, m, u) = -u
_reverse_step(::DiscreteSymmetric, m, step) =
    throw(ArgumentError("move $(move_name(m)) has discrete steps but no reverse_step method"))

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
                    logdet = (vb, x, x′, ctx) -> 0.0, initial_scale = 0.05, invariant = true)

A random-walk move that shifts the `index`-th latent coordinate of the parameter at path
`shift` by the step `u` and replaces the parameter at path `block` with
`compensate(vb, x, x′, ctx)`, where `vb` is its value at `x` and `x′` is `x` with the shift
applied; `compensate` keeps a quantity the likelihood depends on fixed (for example a product
or a pixelwise image). `logdet(vb, x, x′, ctx)` is `log|det ∂x′/∂x|` of the whole map.

Construction checks only that `shift` and `block` exist in `view` and do not overlap.
"""
struct CompensatedMove{C, L} <: AbstractMove
    name::String
    ishift::Int
    block::Tuple
    compensate::C
    logdet::L
    initial_scale::Float64
    invariant::Bool
end

function CompensatedMove(
        name::AbstractString, view::CoordinateView, shift::Tuple, block::Tuple, compensate;
        index::Integer = 1, logdet = (vb, x, x′, ctx) -> 0.0, initial_scale::Real = 0.05,
        invariant::Bool = true
    )
    rs = coords(view, shift)
    1 <= index <= length(rs) || throw(
        ArgumentError("index $index is outside the $(length(rs)) latent coordinates of $shift")
    )
    rb = coords(view, block)
    i = rs[index]
    i in rb && throw(ArgumentError("the shifted coordinate of $shift lies in the block $block"))
    return CompensatedMove(
        String(name), i, block, compensate, logdet, Float64(RandomWalk(initial_scale).initial_scale), invariant
    )
end

step_kind(m::CompensatedMove) = RandomWalk(m.initial_scale)
move_name(m::CompensatedMove) = m.name
is_invariant(m::CompensatedMove) = m.invariant

function propose(m::CompensatedMove, x, u, ctx)
    x′ = _with_coordinate(x, m.ishift, _rget(x, m.ishift) + u)
    vb = value(ctx.view, x, m.block)
    rb = coords(ctx.view, m.block)
    x′[rb] = latent(ctx.view, m.block, m.compensate(vb, x, x′, ctx))
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

_logprior_latent(tbase, post, x) =
    logdensityof(tbase, x) - loglikelihood(post, transform(tbase, x))

"""
    check_move(move, post, θs; space = nothing, rng = Random.default_rng(), τ = nothing,
               rtol = 1e-8, h = 1e-6, logdet_atol = 1e-5, nprior = 0, nsteps = 5,
               ncoords = 5, ks_alpha = 1e-3) -> NamedTuple

Check `move` on the posterior `post` in the base latent space `space` at each constrained
point in `θs`, throwing an error that names the failed property:

  - reversal: proposing with `reverse_step` from the proposal returns to the start, with the
    opposite log-determinant;
  - log-determinant: `logdet` matches `log|det J|` of a central finite-difference Jacobian
    (step `h`) over the coordinates the proposal changes (the map is the identity on the
    others, so they do not contribute);
  - invariance: if `is_invariant(move)`, the log-likelihood does not change (`rtol`).

With `nprior > 0` it also checks that the move kernel leaves the prior invariant: `nprior`
prior draws each take `nsteps` Metropolis–Hastings steps of `move` targeting the prior, and
the result is compared with `nprior` fresh prior draws on up to `ncoords` of the changed
coordinates by a two-sample Kolmogorov–Smirnov test at level `ks_alpha` (Bonferroni over the
coordinates).

The step is drawn with `draw_step(move, rng, τ)`, `τ` defaulting to the move's initial
scale. Returns the largest reversal error, log-determinant error and log-likelihood change,
and the KS statistics and critical value.
"""
function check_move(
        move::AbstractMove, post::VLBIPosterior, θs; space = nothing,
        rng::AbstractRNG = Random.default_rng(), τ = nothing, rtol::Real = 1.0e-8,
        h::Real = 1.0e-6, logdet_atol::Real = 1.0e-5, nprior::Integer = 0, nsteps::Integer = 5,
        ncoords::Integer = 5, ks_alpha::Real = 1.0e-3
    )
    view = CoordinateView(post, space)
    ctx = move_context(view, (move,))
    kind = step_kind(move)
    τ = something(τ, kind isa RandomWalk ? kind.initial_scale : 1.0)
    name = move_name(move)
    tbase = view.tbase
    rev, lderr, dll = 0.0, 0.0, 0.0
    changed = Int[]
    for θ in θs
        x = inverse(tbase, θ)
        step = draw_step(move, rng, τ)
        x′, ld = propose(move, x, step, ctx)
        x″, ld′ = propose(move, x′, reverse_step(move, step), ctx)
        e = maximum(abs, x″ .- x) / (maximum(abs, x) + 1)
        rev = max(rev, e)
        (e <= rtol && isapprox(ld′, -ld; atol = logdet_atol)) || error(
            "move $name is not reversed by reverse_step at step $step: |Δx| = $e, " *
                "logdet $ld then $ld′ (should sum to 0)"
        )
        T = findall(x′ .!= x)
        isempty(T) && error("move $name left the point unchanged at step $step")
        isempty(changed) && (changed = T)
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
        if is_invariant(move)
            l0 = loglikelihood(post, transform(tbase, x))
            l1 = loglikelihood(post, transform(tbase, x′))
            dll = max(dll, abs(l1 - l0))
            abs(l1 - l0) <= rtol * (abs(l0) + 1) || error(
                "move $name changed the log-likelihood from $l0 to $l1; the model is not " *
                    "invariant under it"
            )
        end
    end
    ks, crit = Float64[], NaN
    if nprior > 0
        idx = changed[unique(round.(Int, range(1, length(changed), length = min(ncoords, length(changed)))))]
        crit = sqrt(-log(ks_alpha / length(idx) / 2) / 2) * sqrt(2 / nprior)
        draw() = inverse(tbase, prior_sample(rng, post))
        moved = map(1:nprior) do _
            x = draw()
            ℓ = _logprior_latent(tbase, post, x)
            for _ in 1:nsteps
                x′, ld = propose(move, x, draw_step(move, rng, τ), ctx)
                ℓ′ = _logprior_latent(tbase, post, x′)
                if log(rand(rng)) < ℓ′ - ℓ + ld
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
    return (; reversal = rev, logdet = lderr, loglikelihood = dll, ks, ks_critical = crit)
end
