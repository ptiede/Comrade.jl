# Moves along exact symmetries of instrument models with Gauss–Markov gain chains.

export PhaseSheetMove, flux_gain_move, ChainHyperMove, chain_hyper_moves

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
gain log-amplitudes at path `gains` (a Gauss–Markov chain or per-site constants): a step `u`
on the flux's latent coordinate takes `F → F′` and shifts every gain value by `−c` with
`c = log(F′/F) / power`. When each visibility depends on these parameters only through
`F · exp(power · g)` (e.g. `power = 2` for a common log-amplitude `g` of both stations of a
baseline), the likelihood is unchanged. For a whitened chain with fixed hyperparameters, or
Gaussian per-site constants, the shift is a constant in the latent space, so `logdet = 0`;
`check_move` verifies both.
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

"""
    ChainHyperMove

A [`ComponentwiseRandomWalk`](@ref) move of one fitted hyperparameter field (e.g. `σ` or
`τ`) of a Gauss–Markov chain term, one component per site that fits it; build them with
[`chain_hyper_moves`](@ref). Component `i` steps the latent coordinate of site `i`'s
hyperparameter, `h → h′`, keeps the site's chain values `g` fixed and whitens them again
under `h′`, so the likelihood is unchanged. Its log-determinant is
`log p(g | h′) − log p(g | h) + (‖z′‖² − ‖z‖²)/2` over the site's whitened innovations
`z → z′`, so that with the standard normal latent density its acceptance is that of the
hyperparameter given the chain values. Sites touch disjoint latent coordinates and their
chain densities are separate terms, so the components are independent.
"""
struct ChainHyperMove{P <: Tuple, S} <: AbstractMove
    name::String
    path::P
    sites::Vector{Symbol}
    hcoords::Vector{Int}
    innovations::Vector{Vector{Int}}
    specs::S
    initial_scale::Float64
end

step_kind(m::ChainHyperMove) = ComponentwiseRandomWalk(m.initial_scale, length(m.sites))
move_name(m::ChainHyperMove) = m.name
component_coords(m::ChainHyperMove) = [vcat(h, z) for (h, z) in zip(m.hcoords, m.innovations)]
component_label(m::ChainHyperMove, i) = m.sites[i]

_site_chain_logpdf(specs, g, hp) = sum(spec -> chain_term(spec, g, hp), specs; init = zero(eltype(g)))

function propose(m::ChainHyperMove, x, u::AbstractVector, ctx)
    v = value(ctx.view, x, m.path)
    x′ = copy(x)
    x′[m.hcoords] .+= u
    h′ = value(ctx.view, x′, m.path).hyperparams
    x′[coords(ctx.view, m.path)] = latent(ctx.view, m.path, (params = v.params, hyperparams = h′))
    g = parent(v.params)
    logdet = map(eachindex(m.sites)) do i
        z, z′ = x[m.innovations[i]], x′[m.innovations[i]]
        return _site_chain_logpdf(m.specs[i], g, h′) - _site_chain_logpdf(m.specs[i], g, v.hyperparams) +
            (sum(abs2, z′) - sum(abs2, z)) / 2
    end
    return x′, logdet
end

"""
    chain_hyper_moves(view::CoordinateView, path; initial_scale = 0.1) -> Vector{ChainHyperMove}

One [`ChainHyperMove`](@ref) per fitted hyperparameter field of the Gauss–Markov chain term
at `path` (e.g. `(:instrument, :lg1)`), named `"chain_hyper[term.field]"`, with one
component per site that fits the field.

Errors if `path` is not a chain term with fitted hyperparameters.
"""
function chain_hyper_moves(view::CoordinateView, path::Tuple; initial_scale::Real = 0.1)
    t = node(view, path)
    hasproperty(t, :hnode) || throw(
        ArgumentError("$path is not a Gauss–Markov chain term with fitted hyperparameters")
    )
    r = coords(view, path)
    nh = _node_dimension(t.hnode)
    nh > 0 || throw(ArgumentError("the chain term $path fits no hyperparameters"))
    x0 = randn(Random.Xoshiro(1), dimension(view.tbase)) ./ 3
    v0 = value(view, x0, path)
    h0 = v0.hyperparams
    labels = map(r[1:nh]) do j
        x1 = copy(x0)
        x1[j] += 0.1
        h1 = value(view, x1, path).hyperparams
        moved = [(s, f) for s in keys(h0) for f in keys(h0[s]) if h1[s][f] != h0[s][f]]
        length(moved) == 1 || error(
            "latent coordinate $j of $path moves the hyperparameters $moved; each " *
                "coordinate must move exactly one"
        )
        return only(moved)
    end
    chains = values(t.dists.chains)
    tag = path[end]
    return map(unique(last.(labels))) do f
        idx = findall(l -> last(l) == f, labels)
        sites = first.(labels[idx])
        hcoords = collect(r[idx])
        innovations = map(hcoords) do j
            x1 = copy(x0)
            x1[j] += 0.1
            y = latent(view, path, (params = v0.params, hyperparams = value(view, x1, path).hyperparams))
            k = findall(i -> abs(y[i] - x0[r[i]]) > 1.0e-10 * (abs(x0[r[i]]) + 1), (nh + 1):length(r))
            return collect(r[nh .+ k])
        end
        specs = map(s -> Tuple(c for c in chains if c isa MarkovChainSpec && c.hpsel === Val(s)), sites)
        return ChainHyperMove(
            "chain_hyper[$tag.$f]", path, sites, hcoords, innovations, specs, Float64(initial_scale)
        )
    end
end
