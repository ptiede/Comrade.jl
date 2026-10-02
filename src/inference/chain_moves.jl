# Moves along exact symmetries of instrument models with Gauss–Markov gain chains.

export PhaseSheetMove, flux_gain_move

_is_chain_value(v) = v isa NamedTuple && haskey(v, :params) && haskey(v, :hyperparams)
_chain_params(v) = _is_chain_value(v) ? v.params : v
_with_chain_params(v, p) = _is_chain_value(v) ? (; params = p, hyperparams = v.hyperparams) : p

# Indices of each site's points within a `SiteArray`, in time order.
function _site_points(sa)
    st = sites(sa)
    ts = [t.t0 for t in times(sa)]
    return map(unique(st)) do s
        I = findall(==(s), st)
        s => I[sortperm(ts[I])]
    end
end

# The points of the chain at `path` that the model fixes (reference or initial values): those
# whose value does not move when the latent point does.
function _fixed_chain_points(view::CoordinateView, path)
    n = dimension(view.tbase)
    a = parent(_chain_params(value(view, zeros(n), path)))
    b = parent(_chain_params(value(view, randn(Random.Xoshiro(1), n), path)))
    return a .== b
end

"""
    PhaseSheetMove(post::VLBIPosterior, terms; space = nothing)

Moves between the 2π sheets of the real-line Gauss–Markov phase chains named by `terms`
(instrument parameter names, e.g. `(:gp,)`). The likelihood sees a phase only through
`e^{iφ}`, so adding `±2π` to one site's path from a free point up to the next point the
model fixes (or the end of the path) leaves it unchanged while the chain prior changes.
A step draws a term and site, a free start point and a sign uniformly; the opposite sign
reverses it. The chain's latent block is mapped to its value, shifted and mapped back; with
the hyperparameters unchanged that is a constant latent shift, so `logdet = 0`.

Errors if a term is not an instrument parameter, or has no free point.
"""
struct PhaseSheetMove <: AbstractMove
    paths::Vector{Tuple{Symbol, Symbol}}
    points::Vector{Tuple{Int, Symbol, Vector{Int}, Vector{Int}}}
    fixed::Vector{BitVector}
end

function PhaseSheetMove(post::VLBIPosterior, terms; space = nothing)
    view = CoordinateView(post, space)
    paths = [(:instrument, Symbol(t)) for t in terms]
    isempty(paths) && throw(ArgumentError("PhaseSheetMove needs at least one phase term"))
    v0 = transform(view.tbase, zeros(dimension(view.tbase)))
    points = Tuple{Int, Symbol, Vector{Int}, Vector{Int}}[]
    fixed = BitVector[]
    for (k, path) in enumerate(paths)
        haskey(v0.instrument, path[2]) ||
            throw(ArgumentError("no instrument parameter $(path[2]); found $(keys(v0.instrument))"))
        fx = _fixed_chain_points(view, path)
        push!(fixed, fx)
        for (s, I) in _site_points(_chain_params(v0.instrument[path[2]]))
            free = [j for j in eachindex(I) if !fx[I[j]]]
            isempty(free) || push!(points, (k, s, I, free))
        end
    end
    isempty(points) && throw(ArgumentError("the phase terms $(collect(terms)) have no free point"))
    return PhaseSheetMove(paths, points, fixed)
end

step_kind(::PhaseSheetMove) = DiscreteSymmetric()
move_name(::PhaseSheetMove) = "phase_sheet"
draw_step(m::PhaseSheetMove, rng::AbstractRNG, τ) =
    (p = rand(rng, eachindex(m.points)); (p, rand(rng, m.points[p][4]), rand(rng, (-1, 1))))
reverse_step(::PhaseSheetMove, (p, j, s)) = (p, j, -s)
step_label(m::PhaseSheetMove, (p, _, _)) = (m.paths[m.points[p][1]][2], m.points[p][2])

# The positions of the path `I` a shift from its `j`-th point moves: up to the next fixed one.
function _shift_span(fx, I, j)
    e = findnext(k -> fx[k], I, j + 1)
    return I[j:(isnothing(e) ? lastindex(I) : e - 1)]
end

function propose(m::PhaseSheetMove, x, (p, j, s), ctx)
    k, _, I, _ = m.points[p]
    path = m.paths[k]
    v = value(ctx.view, x, path)
    pa = deepcopy(_chain_params(v))
    parent(pa)[_shift_span(m.fixed[k], I, j)] .+= s * 2π
    x′ = copy(x)
    x′[coords(ctx.view, path)] = latent(ctx.view, path, _with_chain_params(v, pa))
    back = parent(_chain_params(value(ctx.view, x′, path)))
    maximum(abs, back .- parent(pa)) < 1.0e-8 || error(
        "the shifted $(path[2]) path does not round-trip through its transform; the chain " *
            "may be wrapped (a sheet move needs a real-line chain)"
    )
    return x′, zero(eltype(x))
end

"""
    flux_gain_move(view::CoordinateView; flux, gains, power = 2, initial_scale = 0.05)

A [`CompensatedMove`](@ref) trading the scalar flux parameter at path `flux` against the
gain log-amplitude chain at path `gains`: a step `u` on the flux's latent coordinate takes
`F → F′` and shifts every value of the chain by `−c` with `c = log(F′/F) / power`. When
each visibility depends on these parameters only through `F · exp(power · g)` (e.g.
`power = 2` for a common log-amplitude `g` of both stations of a baseline), the likelihood
is unchanged. For a whitened chain with fixed hyperparameters the shift is a constant in
the latent space, so `logdet = 0`; `check_move` verifies both.
"""
function flux_gain_move(
        view::CoordinateView; flux::Tuple, gains::Tuple, power::Real = 2, initial_scale::Real = 0.05
    )
    power > 0 || throw(ArgumentError("power must be positive, got $power"))
    compensate = function (vg, x, x′, ctx)
        c = log(value(ctx.view, x′, flux) / value(ctx.view, x, flux)) / power
        pa = deepcopy(_chain_params(vg))
        parent(pa) .-= c
        return _with_chain_params(vg, pa)
    end
    return CompensatedMove("flux_gain", view, flux, gains, compensate; initial_scale)
end
