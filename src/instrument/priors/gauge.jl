# Finds the parameter entries that must be pinned so that a station phase written as a sum
# of several `ArrayPrior` terms is determined by the data.
#
# A datum is one (site, stamp) pair, where a stamp is one (time, frequency) of the
# observation. Baseline phases depend only on differences of station phases, so the data
# determine the total phase at the sites of a stamp up to one constant `g_τ` per stamp. Each
# datum gives one equation,
#
#     sum_k x_{k,(s,τ)} − g_τ = (data),
#
# where `x_k` runs over the terms of the sum and pinned entries drop out as constants. All
# coefficients are 1, so the system is a hypergraph with one node per unknown (free entries
# and stamp constants) and one hyperedge per datum.

# Union-find over gauge nodes, with path halving and union by rank.
struct GaugeUnionFind
    parent::Vector{Int}
    rank::Vector{Int}
end

GaugeUnionFind(n::Int) = GaugeUnionFind(collect(1:n), zeros(Int, n))

function _find!(uf::GaugeUnionFind, i::Int)
    p = uf.parent
    while p[i] != i
        p[i] = p[p[i]]
        i = p[i]
    end
    return i
end

function _union!(uf::GaugeUnionFind, i::Int, j::Int)
    a, b = _find!(uf, i), _find!(uf, j)
    a == b && return a
    if uf.rank[a] < uf.rank[b]
        a, b = b, a
    end
    uf.parent[b] = a
    uf.rank[a] == uf.rank[b] && (uf.rank[a] += 1)
    return a
end

"""
    default_gauge_preference(entry)

Rank the entries that [`gauge_pins`](@ref) may pin: the entry with the smallest key is the
one it pins, and the largest key is the one it treats as determined when a datum can only
determine one of its entries.

`entry` is a NamedTuple `(; term, index, site, dt, sefd, occupancy)`: the term's name, the
entry's position in the term's `SiteLookup`, its site, the width of its integration time,
the site's SEFD (`Inf` when the array table has no such site), and the number of stamps at
which the site appears.

The key is `(-dt, sefd, -occupancy)`: pin the entry with the longest integration time
first, then the most sensitive site, then the site with the most data.
"""
default_gauge_preference(e) = (-e.dt, e.sefd, -e.occupancy)

# The sorted (time, frequency) stamps of the data and the sites observed at each. Only
# sites connected by baselines share a phase constant, so a (time, frequency) whose
# baselines form several disconnected groups yields one stamp per group, with the same key.
function gauge_stamps(array::AbstractArrayConfiguration)
    T = array[:Ti]
    F = array[:Fr]
    bl = array[:sites]
    keys = Tuple{eltype(T), eltype(F)}[]
    bls_at = Vector{Tuple{Symbol, Symbol}}[]
    index = Dict{Tuple{eltype(T), eltype(F)}, Int}()
    for i in eachindex(T, F, bl)
        key = (T[i], F[i])
        k = get!(index, key) do
            push!(keys, key)
            push!(bls_at, Tuple{Symbol, Symbol}[])
            length(keys)
        end
        push!(bls_at[k], (bl[i][1], bl[i][2]))
    end
    stamps = eltype(keys)[]
    sites_at = Vector{Symbol}[]
    for k in sortperm(keys), group in _connected_sites(bls_at[k])
        push!(stamps, keys[k])
        push!(sites_at, group)
    end
    return stamps, sites_at
end

# The sites of `baselines` grouped into the connected components of the baseline graph,
# each group and the groups themselves in order of first appearance.
function _connected_sites(baselines)
    sites = Symbol[]
    idx = Dict{Symbol, Int}()
    for b in baselines, s in b
        get!(idx, s) do
            push!(sites, s)
            length(sites)
        end
    end
    uf = GaugeUnionFind(length(sites))
    for (s1, s2) in baselines
        _union!(uf, idx[s1], idx[s2])
    end
    groups = Vector{Symbol}[]
    groupof = Dict{Int, Int}()
    for (i, s) in pairs(sites)
        g = get!(groupof, _find!(uf, i)) do
            push!(groups, Symbol[])
            length(groups)
        end
        push!(groups[g], s)
    end
    return groups
end

# For each stamp, the entry of `smap` covering each site there. `searchsorted` over the
# time-sorted stamps is only a prefilter; containment is decided by `in`, whose intervals
# are half open.
function _coverage(smap::SiteLookup, stamps, stamptimes, name)
    cover = [Dict{Symbol, Int}() for _ in eachindex(stamps)]
    for i in eachindex(smap.sites)
        s = smap.sites[i]
        for τ in searchsorted(stamptimes, smap.times[i])
            (stamptimes[τ] ∈ smap.times[i] && last(stamps[τ]) ∈ smap.frequencies[i]) ||
                continue
            prev = get(cover[τ], s, 0)
            prev == 0 || throw(
                ArgumentError(
                    "site $s at $(stamps[τ]) is covered by two entries ($prev and $i) of " *
                        "term $name; its time or frequency segments overlap"
                )
            )
            cover[τ][s] = i
        end
    end
    return cover
end

"""
    gauge_pins(terms, array::AbstractArrayConfiguration; prefer = default_gauge_preference)
    gauge_pins(terms, stamps, sites_at, sefd; prefer = default_gauge_preference)

Return the entries that must be pinned, as a `Vector{Tuple{Symbol, Int}}` of (term name,
entry index), so that a station phase written as the sum of `terms` is fully determined by
the data. An empty result means the parameterization is already identified.

`terms` is a NamedTuple mapping each term's name to `(; smap, fixed)`, where `smap` is the
term's [`SiteLookup`](@ref) and `fixed` holds the entry indices that its referencing scheme
and init pins hold constant. The first form reads the schedule off an array configuration.
The second takes `stamps`, a vector of `(time, frequency)` pairs sorted by time; `sites_at`,
the sites present at each stamp; and `sefd`, a site-keyed NamedTuple or dictionary of system
equivalent flux densities. The sites of one stamp must be connected by baselines: a
`(time, frequency)` whose baselines split into disconnected groups is listed once per group.

Each datum (site `s` at stamp `τ`) gives one equation relating the stamp's phase constant to
the free entries covering `(s, τ)`. An entry covers a datum when its site matches and its
integration time and frequency channel contain the stamp. The solver:

 1. Removes each equation with three or more unknowns that contains an unknown appearing in
    no other equation, and marks that unknown as determined (a stamp constant if there is
    one, otherwise the least preferred entry). This repeats until no equation can be
    removed.
 2. Treats the remaining equations as graph edges: an equation with one unknown determines
    it, and one with two unknowns links them. A remaining equation with three or more
    unknowns throws an `ArgumentError`.
 3. Pins one entry in each connected group of unknowns that no equation determines, chosen
    by `prefer`; see [`default_gauge_preference`](@ref).

A group containing only stamp constants needs no pin.

## Notes
[`InstrumentModel`](@ref) applies this analysis, through its `gaugefix` setting, to the
terms declared with `gauge = :phase`.
"""
function gauge_pins(
        terms, array::AbstractArrayConfiguration; prefer = default_gauge_preference
    )
    stamps, sites_at = gauge_stamps(array)
    tarr = array.tarr
    sefd = NamedTuple{Tuple(tarr.sites)}(Tuple(tarr.SEFD1 .+ tarr.SEFD2))
    return gauge_pins(terms, stamps, sites_at, sefd; prefer)
end

function gauge_pins(terms, stamps, sites_at, sefd; prefer = default_gauge_preference)
    length(stamps) == length(sites_at) || throw(
        DimensionMismatch(
            "each stamp needs its site list: got $(length(stamps)) stamps and " *
                "$(length(sites_at)) site lists"
        )
    )
    names = collect(Symbol, keys(terms))
    smaps = [t.smap for t in values(terms)]
    fixedsets = [Set{Int}(t.fixed) for t in values(terms)]
    isempty(names) && return Tuple{Symbol, Int}[]
    for k in eachindex(smaps, fixedsets)
        all(in(eachindex(smaps[k].sites)), fixedsets[k]) || throw(
            ArgumentError(
                "term $(names[k]) has fixed indices outside its " *
                    "$(length(smaps[k].sites)) entries"
            )
        )
    end

    occupancy = Dict{Symbol, Int}()
    for ss in sites_at, s in ss
        occupancy[s] = get(occupancy, s, 0) + 1
    end

    # Nodes 1:nstamp are the stamp constants, then the entries of each term in turn. Fixed
    # entries get a node that no equation mentions.
    nstamp = length(stamps)
    offsets = Vector{Int}(undef, length(smaps))
    nnodes = nstamp
    for k in eachindex(smaps)
        offsets[k] = nnodes
        nnodes += length(smaps[k].sites)
    end

    stamptimes = first.(stamps)
    issorted(stamptimes) ||
        throw(ArgumentError("the stamps must be sorted by time, then frequency"))
    covers = [_coverage(smaps[k], stamps, stamptimes, names[k]) for k in eachindex(smaps)]

    edges = Vector{Int}[]
    edge_site = Symbol[]
    edge_stamp = Int[]
    for τ in eachindex(stamps, sites_at)
        for s in sites_at[τ]
            e = [τ]
            for k in eachindex(smaps)
                i = get(covers[k][τ], s, 0)
                (i == 0 || i ∈ fixedsets[k]) && continue
                push!(e, offsets[k] + i)
            end
            push!(edges, e)
            push!(edge_site, s)
            push!(edge_stamp, τ)
        end
    end

    incident = [Int[] for _ in 1:nnodes]
    for (ei, e) in pairs(edges), v in e
        push!(incident[v], ei)
    end
    degree = map(length, incident)

    function entry_of(v)
        k = something(findlast(<(v), offsets))
        return k, v - offsets[k]
    end
    function entry_key(v)
        k, i = entry_of(v)
        s = smaps[k].sites[i]
        return prefer(
            (
                term = names[k], index = i, site = s, dt = float(_region(smaps[k].times[i])),
                sefd = float(get(sefd, s, Inf)), occupancy = get(occupancy, s, 0),
            )
        )
    end

    alive = trues(length(edges))
    determined = falses(nnodes)
    queue = [v for v in 1:nnodes if degree[v] == 1]
    while !isempty(queue)
        u = popfirst!(queue)
        (degree[u] == 1 && !determined[u]) || continue
        ei = 0
        for e in incident[u]
            if alive[e]
                ei = e
                break
            end
        end
        # Equations with at most two unknowns are kept for the union-find step.
        (ei == 0 || length(edges[ei]) < 3) && continue
        trigger = _peel_trigger(edges[ei], degree, determined, nstamp, entry_key)
        determined[trigger] = true
        alive[ei] = false
        for v in edges[ei]
            degree[v] -= 1
            (degree[v] == 1 && !determined[v]) && push!(queue, v)
        end
    end

    uf = GaugeUnionFind(nnodes)
    groundnodes = Int[]
    for ei in eachindex(edges)
        alive[ei] || continue
        e = edges[ei]
        if length(e) >= 3
            free = join(("$(names[k])" for k in first.(entry_of.(e[2:end]))), ", ")
            throw(
                ArgumentError(
                    "the datum at site $(edge_site[ei]), stamp $(stamps[edge_stamp[ei]]) " *
                        "leaves $(length(e) - 1) gauge terms free at once ($free). " *
                        "gauge_pins solves parameterizations in which every datum ties one " *
                        "term to the stamp; pin or coarsen one of the terms so that at most " *
                        "one of them is free per datum."
                )
            )
        elseif length(e) == 2
            _union!(uf, e[1], e[2])
        elseif length(e) == 1
            push!(groundnodes, e[1])
        end
    end
    grounded = Set(_find!(uf, v) for v in groundnodes)

    # Undetermined entries grouped by component. Entries that no datum sees are skipped:
    # pinning them identifies nothing.
    components = Dict{Int, Vector{Int}}()
    for k in eachindex(smaps), i in eachindex(smaps[k].sites)
        v = offsets[k] + i
        (i ∈ fixedsets[k] || determined[v] || isempty(incident[v])) && continue
        push!(get!(components, _find!(uf, v), Int[]), v)
    end

    pins = Tuple{Int, Int}[]
    for root in sort!(collect(keys(components)))
        root ∈ grounded && continue
        push!(pins, entry_of(argmin(entry_key, components[root])))
    end
    sort!(pins)
    return [(names[k], i) for (k, i) in pins]
end

# The unknown an equation determines when removed: a stamp constant if one is available,
# otherwise the least preferred entry.
function _peel_trigger(e, degree, determined, nstamp, entry_key)
    best = 0
    for v in e
        (degree[v] == 1 && !determined[v]) || continue
        v <= nstamp && return v
        (best == 0 || entry_key(v) > entry_key(best)) && (best = v)
    end
    return best
end

"""
    fixed_indices(d)

The entry indices an observed instrument prior holds constant: the pins of its referencing
scheme and its init pins.
"""
fixed_indices(d::ObservedArrayPrior) = fixed_indices(d.dists)
fixed_indices(d::ObservedHierarchicalArrayPrior) = fixed_indices(d.dists)
fixed_indices(d::PartiallyConditionedDist) = collect(Int, d.fixed_index)
fixed_indices(d::GaussMarkovChainDist) = collect(Int, d.fixedinds)
fixed_indices(::Dists.Distribution) = Int[]

sitelookup(d::ObservedArrayPrior) = d.sitemap
sitelookup(d::ObservedHierarchicalArrayPrior) = d.sitemap

# The `terms` argument of `gauge_pins` for the observed priors `names`.
function gauge_terms(obs::NamedTuple, names::Tuple)
    return NamedTuple{names}(
        map(names) do n
            d = getproperty(obs, n)
            return (smap = sitelookup(d), fixed = fixed_indices(d))
        end
    )
end

function _pin_description(obs::NamedTuple, name::Symbol, i::Int)
    smap = sitelookup(getproperty(obs, name))
    return "$name[site = $(smap.sites[i]), t = $(_center(smap.times[i]))]"
end

# Check that the reference pins of the `gauge = :phase` terms identify their sum. With
# `gaugefix = :error` a remaining flat direction throws; with `:pin` it is pinned and the
# affected priors are rebuilt.
function fix_phase_gauge(prior::NamedTuple, obs::NamedTuple, array, gaugefix::Symbol)
    names = filter(n -> getproperty(prior, n).gauge === :phase, keys(prior))
    isempty(names) && return obs
    pins = gauge_pins(gauge_terms(obs, names), array)
    isempty(pins) && return obs

    listing = join(("  " * _pin_description(obs, n, i) for (n, i) in pins), "\n")
    gaugefix === :error && throw(
        ArgumentError(
            "the phase parameterization is rank deficient by $(length(pins)): the data " *
                "determine the sum of the gauge = :phase terms only up to that many free " *
                "constants. Pinning\n$listing\nwould remove them. Either set " *
                "gaugefix = :pin on the InstrumentModel to apply these pins " *
                "automatically, or pin the sites yourself with a MultiReference on the " *
                "named term."
        )
    )

    out = obs
    for name in unique(first.(pins))
        inds = [i for (n, i) in pins if n === name]
        p = getproperty(prior, name)
        val = something(reference_value(p.refant), 0.0)
        pinned = ArrayPrior(
            p.default_dist, p.override_dist,
            CompositeReference(p.refant, EntryReference(inds, val)),
            p.phase, p.centroid_station, p.gauge
        )
        out = merge(out, NamedTuple{(name,)}((ObservedArrayPrior(pinned, array),)))
        for i in inds
            @info "gaugefix: pinning $(_pin_description(obs, name, i)) to $val"
        end
    end

    left = gauge_pins(gauge_terms(out, names), array)
    isempty(left) || error(
        "gaugefix = :pin left $(length(left)) phase gauge freedom(s) unfixed: $left"
    )
    return out
end
