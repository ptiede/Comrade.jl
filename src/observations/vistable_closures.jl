export ClosureTable, closure_phases, log_closure_amplitudes

"""
    ClosureTable

Closure phases or log closure amplitudes formed from a [`VisTable`](@ref) of visibilities.
Built by [`closure_phases`](@ref) and [`log_closure_amplitudes`](@ref).

Each closure combines the visibilities of one time stamp and one channel. Its value is
`design * x`, where `x` is the phase (closure phases) or the log amplitude (log closure
amplitudes) of the measured visibilities of the source table, ordered by record and then
channel. The noise covariance is `design * Diagonal((σ ./ |V|).^2) * design'`, the high-SNR
limit, which is block diagonal over time stamps and channels. With this covariance the
likelihood does not depend on which independent set of closures is used.

`measurement(ct)` is a map over `(Node,)` whose coordinates are each closure's `Ti`, `Fr`,
time-stamp index `ts` and `stations`.
"""
struct ClosureTable{K <: Union{ClosurePhases, LogClosureAmplitudes}, M, T <: VisTable, D, C}
    product::K
    measurement::M
    vis::T
    design::D
    covariance::C
end

measurement(ct::ClosureTable) = getfield(ct, :measurement)
domain(ct::ClosureTable) = axisdims(measurement(ct))
info(ct::ClosureTable) = info(getfield(ct, :vis))

"""
    noise(ct::ClosureTable)

The standard deviation of each closure: the square root of the diagonal of its covariance.
"""
noise(ct::ClosureTable) = IntensityMap(sqrt.(diag(getfield(ct, :covariance))), domain(ct))

Base.length(ct::ClosureTable) = length(measurement(ct))
Base.size(ct::ClosureTable) = size(measurement(ct))
Base.size(ct::ClosureTable, d) = size(measurement(ct), d)

"""
    closure_phases(t::VisTable)

An independent set of closure phases of the visibilities in `t`, as a
[`ClosureTable`](@ref). At each time stamp and channel the measured baselines form a graph on
the stations, and its `E - V + C` independent cycles (edges, stations, connected pieces) are
the closure phases. Each comes from one baseline outside the spanning forest of highest
signal-to-noise ratio, closed by the forest path between its stations. A complete array gives
triangles; an incomplete one can also give longer cycles.
"""
closure_phases(t::VisTable) = _closuretable(ClosurePhases(), t)

"""
    log_closure_amplitudes(t::VisTable)

An independent set of log closure amplitudes of the visibilities in `t`, as a
[`ClosureTable`](@ref). At each time stamp and channel, every connected piece of the baseline
graph with `E` baselines and `V` stations gives `E - V` closures, or `E - V + 1` when it has no
odd cycle. Each comes from one baseline outside the spanning forest of highest
signal-to-noise ratio: an even cycle closed by the forest gives alternating signs; an odd one
is joined to a fixed odd cycle of the piece through the forest.
"""
log_closure_amplitudes(t::VisTable) = _closuretable(LogClosureAmplitudes(), t)

_closurevalue(::ClosurePhases, V) = angle(V)
_closurevalue(::LogClosureAmplitudes, V) = log(abs(V))

# One measured baseline of a time stamp and channel: stations `a`, `b`, the column of its
# visibility, and the variance `(σ/|V|)²` of its phase or log amplitude.
struct ClosureEdge
    a::Int
    b::Int
    col::Int
    var::Float64
end

# The spanning forest of smallest variance: parent, edge to the parent, and depth per station.
function _spanningforest(edges)
    root = Dict{Int, Int}()
    find(v) = (r = get!(root, v, v); r == v ? v : (root[v] = find(r)))
    tree = falses(length(edges))
    for k in sortperm(edges; by = e -> e.var)
        ra, rb = find(edges[k].a), find(edges[k].b)
        ra == rb && continue
        root[ra] = rb
        tree[k] = true
    end
    adj = Dict{Int, Vector{Tuple{Int, Int}}}()
    for (k, e) in pairs(edges)
        tree[k] || continue
        push!(get!(adj, e.a, Tuple{Int, Int}[]), (e.b, k))
        push!(get!(adj, e.b, Tuple{Int, Int}[]), (e.a, k))
    end
    parent, pedge, depth, piece = Dict{Int, Int}(), Dict{Int, Int}(), Dict{Int, Int}(), Dict{Int, Int}()
    for v in sort!(unique(vcat(map(e -> e.a, edges), map(e -> e.b, edges))))
        haskey(depth, v) && continue
        depth[v], parent[v], pedge[v], piece[v] = 0, v, 0, v
        queue = [v]
        while !isempty(queue)
            x = popfirst!(queue)
            for (y, k) in get(adj, x, Tuple{Int, Int}[])
                haskey(depth, y) && continue
                depth[y], parent[y], pedge[y], piece[y] = depth[x] + 1, x, k, v
                push!(queue, y)
            end
        end
    end
    return (; tree, parent, pedge, depth, piece)
end

# The forest path from `u` to `v` as steps `(edge, from, to)`.
function _forestpath(f, u, v)
    up, down = Tuple{Int, Int, Int}[], Tuple{Int, Int, Int}[]
    while u != v
        if f.depth[u] >= f.depth[v]
            push!(up, (f.pedge[u], u, f.parent[u]))
            u = f.parent[u]
        else
            push!(down, (f.pedge[v], f.parent[v], v))
            v = f.parent[v]
        end
    end
    return vcat(up, reverse!(down))
end

# The cycle of non-forest edge `k`: along it from `a` to `b`, then back along the forest.
_cycle(f, edges, k) = vcat([(k, edges[k].a, edges[k].b)], _forestpath(f, edges[k].b, edges[k].a))

function _addterm!(c, col, x)
    c[col] = get(c, col, 0) + x
    return c
end

# Phases along the steps, oriented by each step's direction.
function _phaseterms(edges, steps)
    c = Dict{Int, Int}()
    for (k, from, _) in steps
        _addterm!(c, edges[k].col, edges[k].a == from ? 1 : -1)
    end
    return c
end

# Log amplitudes along the steps with alternating signs, scaled by `w`.
function _alternating!(c, edges, steps, w)
    for (i, (k, _, _)) in pairs(steps)
        _addterm!(c, edges[k].col, w * (isodd(i) ? 1 : -1))
    end
    return c
end

function _closures(::ClosurePhases, edges)
    f = _spanningforest(edges)
    return [_phaseterms(edges, _cycle(f, edges, k)) for k in eachindex(edges) if !f.tree[k]]
end

function _closures(::LogClosureAmplitudes, edges)
    f = _spanningforest(edges)
    odd = Dict{Int, Int}()
    for k in eachindex(edges)
        e = edges[k]
        !f.tree[k] && iseven(f.depth[e.a] + f.depth[e.b]) && get!(odd, f.piece[e.a], k)
    end
    out = Dict{Int, Int}[]
    for k in eachindex(edges)
        f.tree[k] && continue
        e = edges[k]
        c = _alternating!(Dict{Int, Int}(), edges, _cycle(f, edges, k), 1)
        if iseven(f.depth[e.a] + f.depth[e.b])
            k == odd[f.piece[e.a]] && continue
            o = odd[f.piece[e.a]]
            path = _forestpath(f, e.a, edges[o].a)
            _alternating!(c, edges, _cycle(f, edges, o), isodd(length(path)) ? 1 : -1)
            _alternating!(c, edges, path, -2)
        end
        push!(out, c)
    end
    return out
end

# Every closure must cancel the station terms: signed sums for phases, unsigned for amplitudes.
function _checkclosure(product, edges, c)
    bycol = Dict(e.col => e for e in edges)
    total = Dict{Int, Int}()
    for (col, x) in c
        e = bycol[col]
        _addterm!(total, e.a, x)
        _addterm!(total, e.b, product isa ClosurePhases ? -x : x)
    end
    all(iszero, values(total)) || error("closure $c does not cancel the station terms")
    return nothing
end

function _closuretable(product, t::VisTable)
    eltype(measurement(t)) <: Complex ||
        throw(ArgumentError("closures are formed from a table of scalar complex visibilities, got one of $(eltype(measurement(t)))"))
    d = domain(t)
    measured = _measured(t)
    keep = sort!(findall(measured); by = I -> (I[1], I[2]))
    column = zeros(Int, size(measured))
    column[keep] .= eachindex(keep)
    V, σ = _dense(measurement(t)), _dense(noise(t))
    names = info(t).stations.name

    rows, cols, vals = Int[], Int[], Float64[]
    points = (; Ti = Float64[], Fr = Float64[], ts = Int[], stations = Vector{Symbol}[])
    for ts in unique(d.ts), ch in axes(measured, 2)
        nodes = findall(n -> d.ts[n] == ts && measured[n, ch], eachindex(d.ts))
        seen = Set{Tuple{Int, Int}}()
        edges = map(nodes) do n
            a, b = d.ant1[n], d.ant2[n]
            minmax(a, b) in seen &&
                throw(ArgumentError("baseline $(names[a])-$(names[b]) appears twice at time $(d.Ti[n]) and channel $ch"))
            push!(seen, minmax(a, b))
            ClosureEdge(a, b, column[n, ch], abs2(σ[n, ch] / abs(V[n, ch])))
        end
        isempty(edges) && continue
        for c in _closures(product, edges)
            _checkclosure(product, edges, c)
            r = length(points.Ti) + 1
            for (col, x) in c
                iszero(x) && continue
                push!(rows, r); push!(cols, col); push!(vals, x)
            end
            ids = sort!(unique(vcat([e.a for e in edges if haskey(c, e.col) && !iszero(c[e.col])], [e.b for e in edges if haskey(c, e.col) && !iszero(c[e.col])])))
            push!(points.Ti, d.Ti[first(nodes)]); push!(points.Fr, info(t).channels.frequency[ch])
            push!(points.ts, ts); push!(points.stations, names[ids])
        end
    end
    isempty(points.Ti) && throw(ArgumentError("no closures can be formed from this table"))

    design = sparse(rows, cols, vals, length(points.Ti), length(keep))
    x = _closurevalue.(Ref(product), V[keep])
    cov = design * Diagonal(abs2.(σ[keep] ./ abs.(V[keep]))) * design'
    dom = StructuredDomain(
        (Node(length(points.Ti)),); header = info(t),
        Ti = points.Ti => (:Node,), Fr = points.Fr => (:Node,), ts = points.ts => (:Node,),
        stations = points.stations => (:Node,),
    )
    return ClosureTable(product, IntensityMap(design * x, dom), t, design, cov)
end

function Base.show(io::IO, ::MIME"text/plain", ct::ClosureTable{K}) where {K}
    a = info(ct)
    println(io, "ClosureTable{", nameof(K), "}")
    println(io, "  source:   ", a.source.name, " (", a.epoch, ")")
    return print(io, "  closures: ", length(ct))
end
