# Station cells: one per station, time stamp and channel of a `VisTable`. Instrument terms are
# evaluated once per cell and gathered by each visibility's two stations.

"""
    CellLayout

The station cells of a [`VisTable`](@ref) and the index maps the instrument model reads:

- `coords`: per cell, the station index `ant`, `Ti`, `Fr`, `elevation`, `parallactic`,
  `azimuth`, and the direction `ra`, `dec` of the observed source;
- `cell1`, `cell2`: per value `(record, channel)`, the cells of its first and second station;
- `entries`: per instrument parameter, the entry of its site array at each cell;
- `response`: per cell, the matrix `U` from the global sky basis to the station's local ideal
  feed basis, as its four column-major component vectors, or `nothing` when the data need
  none (Stokes I).
"""
struct CellLayout{C, I, E, R}
    coords::C
    cell1::I
    cell2::I
    entries::E
    response::R
end

EnzymeRules.inactive_type(::Type{<:CellLayout}) = true

# The cells `(station, time stamp, channel)` of every record's two stations, sorted, and the
# cells of each value.
function _cells(t::VisTable)
    d = domain(t)
    chs = axes(measurement(t), 2)
    key1 = [(d.ant1[n], d.ts[n], ch) for n in eachindex(d.ts), ch in chs]
    key2 = [(d.ant2[n], d.ts[n], ch) for n in eachindex(d.ts), ch in chs]
    cellkeys = sort!(unique!(vcat(vec(key1), vec(key2))))
    index = Dict(k => i for (i, k) in pairs(cellkeys))
    return cellkeys, getindex.(Ref(index), key1), getindex.(Ref(index), key2)
end

function _cellcoords(t::VisTable, keys, shift)
    τ = times(t) .+ shift
    freqs = info(t).channels.frequency
    g = stations(t)
    el, par, az = parent(g.elevation), parent(g.parallactic), parent(g.azimuth)
    return StructVector(
        (;
            ant = [k[1] for k in keys], Ti = [τ[k[2]] for k in keys], Fr = [freqs[k[3]] for k in keys],
            elevation = [el[k[1], k[2]] for k in keys], parallactic = [par[k[1], k[2]] for k in keys],
            azimuth = [az[k[1], k[2]] for k in keys],
            ra = fill(deg2rad(info(t).source.ra), length(keys)), dec = fill(deg2rad(info(t).source.dec), length(keys)),
        )
    )
end

# The entry of the site array `x` of the parameter `p` holding each cell: the entry of the
# cell's station whose time segment holds the cell's time and whose frequency segment holds its
# frequency. Within one station and frequency segment the time segments do not overlap, so one
# sorted search finds it. A cell with no measured value, or of a station not in `readers` (the
# stations whose term reads `p`), gets the first entry, whose value it never uses.
function _cellentries(p, x::SiteArray, names, coords, measured, readers)
    st, tm, fr = sites(x), times(x), frequencies(x)
    start(i) = tm[i].t0 - tm[i].dt / 2
    groups = Dict{Any, Vector{Int}}()
    for i in eachindex(st)
        push!(get!(groups, (st[i], fr[i]), Int[]), i)
    end
    bysite = Dict{Any, Vector{Tuple{Any, Vector{Int}, Vector{Float64}}}}()
    for ((s, f), e) in groups
        sort!(e; by = start)
        push!(get!(bysite, s, Tuple{Any, Vector{Int}, Vector{Float64}}[]), (f, e, map(start, e)))
    end
    return map(eachindex(coords.ant)) do c
        s = names[coords.ant[c]]
        t, f = coords.Ti[c], coords.Fr[c]
        found = 0
        for (seg, e, starts) in get(bysite, s, ())
            f in seg || continue
            k = searchsortedlast(starts, t)
            (k >= firstindex(e) && t in tm[e[k]]) || continue
            found == 0 || throw(ArgumentError("station $s at time $t h and frequency $f Hz lies in two segments of the instrument parameter $p"))
            found = e[k]
        end
        if found == 0
            measured[c] && s in readers &&
                throw(ArgumentError("the instrument parameter $p has no entry for station $s at time $t h and frequency $f Hz"))
            return firstindex(st)
        end
        return found
    end
end

# `U` from the global basis `refbasis` to the local ideal feed basis of a station whose two feeds
# have the nominal polarizations `feeds` and are turned against the sky by the angles `φ`: row
# `i` is that of the station's basis rotated by `φ[i]`.
function _response(refbasis, feeds::NTuple{2, Any}, φ::NTuple{2, Real})
    basis = PolBasis{typeof(feeds[1]), typeof(feeds[2])}()
    T = basis_transform(refbasis, basis)
    row(i) = (rotation(φ[i], basis) * T)[i, :]
    r1, r2 = row(1), row(2)
    return SMatrix{2, 2, ComplexF64}(r1[1], r2[1], r1[2], r2[2])
end

# The response of every cell, from its station's mount, feeds and receptor angles.
function _responses(t::VisTable, coords, refbasis)
    st = info(t).stations
    return map(coords.ant, coords.elevation, coords.parallactic, coords.azimuth) do a, el, par, az
        feeds = st.feeds[a]
        length(feeds) == 2 || throw(
            ArgumentError("station $(st.name[a]) has $(length(feeds)) feeds; the instrument model supports stations with two")
        )
        any(isnan, st.feed_angles[a]) &&
            throw(ArgumentError("station $(st.name[a]) has unknown receptor angles; set its `feed_angles`"))
        rot = feed_rotation(st.mount[a], par, el, az)
        return _response(refbasis, feeds, map(ψ -> rot + ψ, st.feed_angles[a]))
    end
end

# `shift`: the hours from the epoch of `t` to that of the instrument's site arrays.
function CellLayout(t::VisTable, params::NamedTuple, refbasis, readers::NamedTuple; shift = 0.0)
    cellkeys, cell1, cell2 = _cells(t)
    coords = _cellcoords(t, cellkeys, shift)
    measured = falses(length(cellkeys))
    meas = _measured(t)
    for I in CartesianIndices(meas)
        meas[I] || continue
        measured[cell1[I]] = true
        measured[cell2[I]] = true
    end
    names = info(t).stations.name
    entries = NamedTuple{keys(params)}(
        map((p, x, r) -> _cellentries(p, x, names, coords, measured, r), keys(params), values(params), values(readers[keys(params)]))
    )
    response = _needsresponse(t.pol) ? _matrixcolumns(_responses(t, coords, refbasis)) : nothing
    return CellLayout(coords, cell1, cell2, entries, response)
end

_needsresponse(::IPol) = false
_needsresponse(_) = true

# A vector of 2×2 matrices as their four column-major component vectors.
_matrixcolumns(m::AbstractVector{<:StaticArraysCore.StaticMatrix}) = ntuple(k -> getindex.(m, k), Val(4))
