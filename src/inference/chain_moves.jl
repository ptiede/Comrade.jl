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

A move that adds `±2π` to one site's real-line Gauss–Markov phase chain, from a free point up
to the next point the model fixes (or the end of the chain). `terms` names the instrument
phase terms, e.g. `(:gp,)`. The likelihood depends on a phase only through `e^{iφ}`, so it is
unchanged; the chain prior changes. A step picks a term, site, starting point and sign
uniformly, and the opposite sign reverses it. The shift is a constant in latent space, so the
log-determinant is 0.

Throws if a term is not an instrument parameter, has no free point, or is a wrapped chain.
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
    m = PhaseSheetMove(paths, points, fixed)
    _check_real_line(m, view)
    return m
end

# Throws unless a 2π shift of each site's path survives the round trip through the chain's
# transform, which fails for a wrapped chain.
function _check_real_line(m::PhaseSheetMove, view::CoordinateView)
    x = zeros(dimension(view.tbase))
    ctx = move_context(view, (m,))
    for (p, (k, s, I, free)) in pairs(m.points)
        path = m.paths[k]
        x′, _ = propose(m, x, (p, first(free), 1), ctx)
        shifted = copy(parent(_chain_params(value(view, x, path))))
        shifted[_shift_span(m.fixed[k], I, first(free))] .+= 2π
        back = parent(_chain_params(value(view, x′, path)))
        maximum(abs, back .- shifted) < sqrt(eps()) * (maximum(abs, shifted) + 1) || throw(
            ArgumentError(
                "a 2π shift of site $s of $(path[2]) does not round-trip through its " *
                    "transform; the chain may be wrapped (a sheet move needs a real-line chain)"
            )
        )
    end
    return nothing
end

function step_distribution(m::PhaseSheetMove)
    draw = function (rng)
        p = rand(rng, eachindex(m.points))
        return p, rand(rng, m.points[p][4]), rand(rng, (-1, 1))
    end
    return SymmetricDiscreteStep(draw, ((p, j, s),) -> (p, j, -s))
end
move_name(::PhaseSheetMove) = "phase_sheet"
likelihood_invariant(::PhaseSheetMove) = true
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
    return x′, zero(eltype(x))
end

"""
    flux_gain_move(view::CoordinateView; flux, gains, power = 2)

A [`CompensatedMove`](@ref) that trades the scalar flux at path `flux` against the gain
log-amplitudes at path `gains` (a Gauss–Markov chain or per-site constants). A step changes
the flux from `F` to `F′` and subtracts `log(F′/F) / power` from every gain. The likelihood is
unchanged when each visibility depends on these only through `F · exp(power · g)`, e.g.
`power = 2` when both stations of a baseline share the log-amplitude `g`. The
log-determinant is 0, which holds when the gain shift is a constant in latent space: a
whitened chain with fixed hyperparameters, or Gaussian per-site constants.
"""
function flux_gain_move(
        view::CoordinateView; flux::Tuple, gains::Tuple, power::Real = 2
    )
    power > 0 || throw(ArgumentError("power must be positive, got $power"))
    compensate = function (vg, x, x′, ctx)
        c = log(value(ctx.view, x′, flux) / value(ctx.view, x, flux)) / power
        pa = deepcopy(_chain_params(vg))
        parent(pa) .-= c
        return _with_chain_params(vg, pa)
    end
    return CompensatedMove("flux_gain", view, flux, gains, compensate; likelihood_invariant = true)
end

"""
    ChainHyperMove

A move of one fitted hyperparameter field (e.g. `σ` or `τ`) of a Gauss–Markov chain term, with
one [`ComponentwiseNormalStep`](@ref) component per site; build it with
[`chain_hyper_moves`](@ref). Component `i` changes the hyperparameter `h` of site `i` and
recomputes the site's whitened chain coordinates `z`, so the chain values `g`, and with them
the likelihood, stay the same. The log-determinant
`log p(g | h′) − log p(g | h) + (‖z′‖² − ‖z‖²)/2` makes the acceptance that of the
hyperparameter given the chain values.
"""
struct ChainHyperMove{P <: Tuple, S} <: AbstractMove
    name::String
    path::P
    sites::Vector{Symbol}
    hcoords::Vector{Int}
    innovations::Vector{Vector{Int}}
    specs::S
end

step_distribution(m::ChainHyperMove) = ComponentwiseNormalStep(length(m.sites))
move_name(m::ChainHyperMove) = m.name
likelihood_invariant(::ChainHyperMove) = true
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

# The latent coordinates of each chain spec: after the hyperparameters, one contiguous
# block per spec in spec order.
function _chain_blocks(t, chains, r)
    dims = [_chain_block_dim(t, c) for c in chains]
    sum(dims) == length(r) || error(
        "the chain specs take $(sum(dims)) latent coordinates, but the chain term has $(length(r))"
    )
    starts = cumsum(dims) .- dims
    return [r[s .+ (1:d)] for (s, d) in zip(starts, dims)]
end

_chain_block_dim(t::StdWhitenedHierarchicalTransform, c) = _std_dim(c, t.space)
_chain_block_dim(::WhitenedHierarchicalTransform, c) = _flat_dim(c)

"""
    chain_hyper_moves(view::CoordinateView, path) -> Vector{ChainHyperMove}

One [`ChainHyperMove`](@ref) per fitted hyperparameter field of the Gauss–Markov chain term at
`path` (e.g. `(:instrument, :lg1)`), named `"chain_hyper[term.field]"`.

Throws if `path` is not a chain term with fitted hyperparameters, or if a site's chain is
grouped or wrapped non-centered: the move needs the chain density of each site.
"""
function chain_hyper_moves(view::CoordinateView, path::Tuple)
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
    blocks = _chain_blocks(t, chains, r[(nh + 1):end])
    tag = path[end]
    return map(unique(last.(labels))) do f
        idx = findall(l -> last(l) == f, labels)
        sites = first.(labels[idx])
        hcoords = collect(r[idx])
        ix = map(s -> [i for (i, c) in pairs(chains) if c isa MarkovChains && c.sites === Val((s,))], sites)
        specs = map(i -> chains[i], ix)
        for (s, sp) in zip(sites, specs)
            isempty(sp) && throw(
                ArgumentError(
                    "site $s of $path has no per-site chain whose density chain_hyper_moves can " *
                        "evaluate (a grouped chain?); the move needs log p(g | h) of the site"
                )
            )
            any(_fused_flat, sp) && throw(
                ArgumentError(
                    "site $s of $path is a wrapped non-centered chain, whose chain_term is not " *
                        "its full density; the move needs log p(g | h) of the site"
                )
            )
        end
        innovations = map(i -> reduce(vcat, blocks[i]; init = Int[]), ix)
        return ChainHyperMove("chain_hyper[$tag.$f]", path, sites, hcoords, innovations, specs)
    end
end
