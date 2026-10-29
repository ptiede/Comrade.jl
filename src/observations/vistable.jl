export ArrayInfo, Source, VisTable, Correlation, Correlations, IPol, amplitudes, weight, info, stations, baselines, drop,
    measurement, noise, domain, sites, beamsize, flag, add_fractional_noise, reset_mounts!

using PolarizedTypes: IPol

"""
    Source(name, ra, dec)

The observed source: its `name`, a `Symbol`, and its J2000 right ascension `ra` and
declination `dec` in degrees.
"""
struct Source
    name::Symbol
    ra::Float64
    dec::Float64
    Source(name, ra, dec) = new(Symbol(name), ra, dec)
end

"""
    ArrayInfo(; stations, channels, scans, source, epoch, tracks = nothing)

Metadata of one observation: the stations, spectral channels, scans, tracks, source and epoch.
It is the header of a [`VisTable`](@ref)'s domain and stays on the host.

- `stations`: a table (any Tables.jl source) with one row per station and the columns
  - `name`: a `Symbol`;
  - `position`: the ITRF position in meters, three numbers;
  - `mount`: the mount, an `AbstractMount` such as `MountAltAz()` or `MountNasmythR()`, whose
    [`feed_rotation`](@ref) turns the feeds against the sky;
  - `feeds`: the nominal polarization of each feed in the station's own frame, e.g.
    `(RPol(), LPol())` or `(XPol(), YPol())`;
  - `feed_angles`: the fixed angle in radians by which each feed is turned from its nominal
    orientation on the dish, one per feed; ideal orthogonal feeds share one angle. The mount
    adds the rotation that changes as the station tracks.
  Any other column is kept as given.
- `channels`: a table with the columns `frequency` and `width`, in Hz, sorted by frequency.
- `scans`: a table with the columns `start` and `stop`, in hours since `epoch`; each scan is
  the closed interval `[start, stop]`, with `start < stop`.
- `tracks`: a table with the columns `start` and `stop`, in hours since `epoch`: the observing
  sessions, e.g. one per night, in order and disjoint, each holding whole scans. `nothing`, the
  default, is one track from the first scan's start to the last scan's stop.
- `source`: the [`Source`](@ref).
- `epoch`: the UTC `Date` whose midnight is time zero.
"""
struct ArrayInfo{S, C, K, T} <: ComradeBase.AbstractHeader{Float64, NamedTuple{(), Tuple{}}}
    stations::S
    channels::C
    scans::K
    tracks::T
    source::Source
    epoch::Date
end

function ArrayInfo(; stations, channels, scans, source, epoch, tracks = nothing)
    cols = _requirecolumns(stations, (:name, :position, :mount, :feeds, :feed_angles), "stations")
    position = map(p -> SVector{3, Float64}(p...), cols.position)
    feeds = map(Tuple, cols.feeds)
    feed_angles = map(a -> Float64.(Tuple(a)), cols.feed_angles)
    st = StructVector(merge(cols, (; name = Symbol.(cols.name), position, feeds, feed_angles)))
    allunique(st.name) || throw(ArgumentError("station names must be unique, got $(st.name)"))
    for i in eachindex(st.name)
        f, a = st.feeds[i], st.feed_angles[i]
        all(p -> p isa PolarizedTypes.ElectricFieldBasis, f) ||
            throw(ArgumentError("station $(st.name[i]): feeds are polarizations like RPol() or XPol(), got $f"))
        length(f) == length(a) ||
            throw(DimensionMismatch("station $(st.name[i]) has $(length(f)) feeds and $(length(a)) feed angles"))
        st.mount[i] isa AbstractMount ||
            throw(ArgumentError("station $(st.name[i]): the mount must be an AbstractMount, got $(st.mount[i])"))
    end

    ch = StructVector(_requirecolumns(channels, (:frequency, :width), "channels"))
    isempty(ch) && throw(ArgumentError("an observation needs at least one channel"))
    issorted(ch.frequency) && allunique(ch.frequency) ||
        throw(ArgumentError("channel frequencies must be strictly increasing, got $(ch.frequency)"))

    sc = StructVector(_requirecolumns(scans, (:start, :stop), "scans"))
    for k in eachindex(sc.start, sc.stop)
        a, b = sc.start[k], sc.stop[k]
        a < b || throw(ArgumentError("scan $k is [$a, $b] h; a scan spans at least one integration, so start < stop"))
    end
    tr = StructVector(_requirecolumns(something(tracks, (; start = [minimum(sc.start)], stop = [maximum(sc.stop)])), (:start, :stop), "tracks"))
    _checktracks(tr, sc)

    source isa Source || throw(ArgumentError("source must be a Source(name, ra, dec), got $(repr(source))"))
    return ArrayInfo(st, ch, sc, tr, source, Date(epoch))
end

function _checktracks(tr, sc)
    isempty(tr) && throw(ArgumentError("an observation needs at least one track"))
    for k in eachindex(tr.start, tr.stop)
        tr.start[k] < tr.stop[k] || throw(ArgumentError("track $k is [$(tr.start[k]), $(tr.stop[k])] h; a track needs start < stop"))
        k == firstindex(tr.start) || tr.stop[k - 1] < tr.start[k] ||
            throw(ArgumentError("the tracks must be in order and disjoint; track $(k - 1) ends at $(tr.stop[k - 1]) h and track $k starts at $(tr.start[k]) h"))
    end
    for k in eachindex(sc.start, sc.stop)
        any(i -> tr.start[i] <= sc.start[k] && sc.stop[k] <= tr.stop[i], eachindex(tr.start, tr.stop)) ||
            throw(ArgumentError("scan $k, [$(sc.start[k]), $(sc.stop[k])] h, is not inside a track"))
    end
    return nothing
end

function _requirecolumns(table, required, what)
    cols = Tables.columntable(table)
    missingcols = filter(c -> !haskey(cols, c), required)
    isempty(missingcols) ||
        throw(ArgumentError("$what needs the columns $(required); missing $(missingcols)"))
    return cols
end

"""
    stationindex(info::ArrayInfo, s)

The position of station `s` (a `Symbol` name or an integer index) in `info.stations`.
"""
function stationindex(info::ArrayInfo, s::Symbol)
    i = findfirst(==(s), info.stations.name)
    isnothing(i) && throw(ArgumentError("station $s is not in the array; stations are $(info.stations.name)"))
    return i
end
function stationindex(info::ArrayInfo, s::Integer)
    s in eachindex(info.stations) ||
        throw(ArgumentError("station index $s is out of range $(eachindex(info.stations))"))
    return Int(s)
end

function Base.show(io::IO, info::ArrayInfo)
    return print(
        io, "ArrayInfo(", info.source.name, ", ", info.epoch, ", ", length(info.stations),
        " stations, ", length(info.channels), " channels, ", length(info.scans), " scans, ",
        length(info.tracks), " track", length(info.tracks) == 1 ? "" : "s", ")"
    )
end

"""
    Correlation(i, j)

One correlation product of a baseline: feed `i` of the first station correlated with feed `j`
of the second. A station has two feeds, `1` and `2`, whose nominal polarizations are its
`feeds` in [`ArrayInfo`](@ref), so `Correlation(1, 1)` of two stations with
`(RPol(), LPol())` feeds is RR. Real feeds are not pure, which the instrument model accounts
for.
"""
struct Correlation{I, J}
    function Correlation{I, J}() where {I, J}
        I in (1, 2) && J in (1, 2) ||
            throw(ArgumentError("a station has the feeds 1 and 2; got the correlation ($I, $J)"))
        return new{I, J}()
    end
end
Correlation(i::Integer, j::Integer) = Correlation{Int(i), Int(j)}()

"""
    Correlations()

All four correlation products of a baseline as a 2×2 matrix, whose entry `(i, j)` is
[`Correlation(i, j)`](@ref Correlation).
"""
struct Correlations end

# What the values of a `VisTable` are.
const TablePolarization = Union{IPol, Correlation, Correlations}

# The shape of one value.
_valueshape(::Union{IPol, Correlation}) = ()
_valueshape(::Correlations) = (2, 2)

_polname(::IPol) = "Stokes I"
_polname(::Correlation{I, J}) where {I, J} = "correlation ($I, $J)"
_polname(::Correlations) = "all correlations"

"""
    VisTable(pol, V, σ; u, v, Ti, ant1, ant2, info, weight = nothing,
             elevation = nothing, parallactic = nothing)

A table of complex visibilities `V` and their weights `w = 1/σ²` over one visibility domain with
the dims `(Node, Fr)`. A `Node` is one record, a baseline at one time; `Fr` holds the channels
of `info` (an [`ArrayInfo`](@ref)).

`pol` states what each value is: `IPol()` for Stokes I, the mean of the parallel-hand products;
a [`Correlation(i, j)`](@ref Correlation) for one product of the two stations' feeds; or
[`Correlations()`](@ref Correlations) for all four products as a 2×2 matrix. With `n` records
and `nch` channels, `V` and `σ` have the size `(n, nch)`, or `(n, nch, 2, 2)` for
`Correlations()` with the feeds of the first and second station in the last two dims. With one
channel, `(n,)` and `(n, 2, 2)` work too. Instead of `σ`, the weights can be given as
`weight = w`. Amplitudes are taken from a table with [`amplitudes`](@ref).

Per record: `u`, `v`, the baseline in meters; `Ti`, the time in hours since `info.epoch`;
`ant1`, `ant2`, the two stations, as names or indices into `info.stations`.

A value without a measurement gets weight `0`: a non-finite or non-positive `σ`, a non-finite
`V`, or a hand missing from a coherency. Its measurement is stored as `0`.

The records are sorted by time. The source elevation and parallactic angle of every station at
every time stamp are computed from the station positions, or taken from the matrices
`elevation` and `parallactic` (stations × time stamps) when given.

See [`measurement`](@ref), [`weight`](@ref), [`noise`](@ref), [`domain`](@ref),
[`info`](@ref), [`stations`](@ref) and [`baselines`](@ref).
"""
struct VisTable{P, M, W, G}
    pol::P
    measurement::M
    weight::W
    stations::G
end

function VisTable(
        pol::TablePolarization, V, σ = nothing;
        u, v, Ti, ant1, ant2, info::ArrayInfo, weight = nothing,
        elevation = nothing, parallactic = nothing
    )
    isnothing(σ) == isnothing(weight) &&
        throw(ArgumentError("give either the noise σ or `weight`, not both or neither"))
    Base.require_one_based_indexing(V, u, v, Ti, ant1, ant2)
    n = length(u)
    all(==(n), map(length, (v, Ti, ant1, ant2))) ||
        throw(DimensionMismatch("u, v, Ti, ant1 and ant2 need one entry per record; got lengths $(map(length, (u, v, Ti, ant1, ant2)))"))
    nch = length(info.channels)

    Vd = _datashape(pol, V, n, nch, "V")
    w = isnothing(weight) ? _weights(Vd, _datashape(pol, σ, n, nch, "σ")) :
        _checkweights(Vd, _datashape(pol, weight, n, nch, "weight"))
    Vd = ifelse.(w .> 0, Vd, zero(eltype(Vd)))

    a1 = map(s -> stationindex(info, s), ant1)
    a2 = map(s -> stationindex(info, s), ant2)
    any(a1 .== a2) && throw(ArgumentError("record $(findfirst(a1 .== a2)) joins a station to itself"))

    perm = sortperm(collect(zip(Ti, a1, a2)))
    times = sort!(unique(Ti))
    t = Float64.(Ti[perm])
    ts = searchsortedfirst.(Ref(times), t)

    dom = StructuredDomain(
        (Node(n), Fr(info.channels.frequency)); header = info,
        u = Float64.(u[perm]) => (:Node,), v = Float64.(v[perm]) => (:Node,), Ti = t => (:Node,),
        ts = ts => (:Node,), ant1 = a1[perm] => (:Node,), ant2 = a2[perm] => (:Node,),
    )
    meas = IntensityMap(_elements(pol, selectdim(Vd, 1, perm)), dom)
    wmap = IntensityMap(_elements(pol, selectdim(w, 1, perm)), dom)
    return VisTable(pol, meas, wmap, _stationgrid(info, times, elevation, parallactic))
end

function _datashape(pol, x, n, nch, what)
    fd = _valueshape(pol)
    if size(x) == (n, nch, fd...)
        y = x
    elseif nch == 1 && size(x) == (n, fd...)
        y = reshape(x, n, 1, fd...)
    else
        throw(DimensionMismatch("$what has size $(size(x)); expected $((n, nch, fd...)) for $n records and $nch channels"))
    end
    y = coalesce.(y, NaN)
    return what == "V" ? complex(float(y)) : float(y)
end

_weights(V, σ) = @. ifelse(isfinite(V) & isfinite(σ) & (σ > 0), inv(σ^2), zero(σ))

function _checkweights(V, w)
    all(x -> isfinite(x) && x >= 0, w) || throw(ArgumentError("weights must be finite and non-negative"))
    return @. ifelse(isfinite(V), w, zero(w))
end

_elements(::Correlations, x) = FieldDimArray{SMatrix{2, 2}}(collect(x))
_elements(_, x) = collect(x)

function _stationgrid(info::ArrayInfo, times, elevation, parallactic)
    isnothing(elevation) == isnothing(parallactic) ||
        throw(ArgumentError("give both `elevation` and `parallactic`, or neither"))
    sz = (length(info.stations), length(times))
    el, par, az = _geometry(info, times)
    if !isnothing(elevation)
        size(elevation) == sz && size(parallactic) == sz ||
            throw(DimensionMismatch("elevation and parallactic need the size (stations, time stamps) = $sz; got $(size(elevation)) and $(size(parallactic))"))
        el, par = Float64.(elevation), Float64.(parallactic)
    end
    gd = (Ant(info.stations.name), Ts(times))
    return DD.DimStack((; elevation = DD.DimArray(el, gd), parallactic = DD.DimArray(par, gd), azimuth = DD.DimArray(az, gd)))
end

function _geometry(info::ArrayInfo, times)
    jd0 = datetime2julian(DateTime(info.epoch))
    st = info.stations
    el, par, az = (zeros(length(st), length(times)) for _ in 1:3)
    for (j, t) in pairs(times), i in eachindex(st.position)
        el[i, j], par[i, j], az[i, j] = horizon_angles(st.position[i], info.source.ra, info.source.dec, jd0 + t / 24)
    end
    return el, par, az
end

"""
    measurement(t::VisTable)

The measurements of `t`: a map over `(Node, Fr)`, the records and the channels, whose
elements are complex visibilities, amplitudes or 2×2 coherency matrices. Unmeasured values
are `0`.
"""
measurement(t::VisTable) = getfield(t, :measurement)

"""
    weight(t::VisTable)

The weights `1/σ²` of the measurements of `t`, `0` where nothing was measured.
"""
weight(t::VisTable) = getfield(t, :weight)

"""
    noise(t::VisTable)

The noise `σ = 1/√w` of the measurements of `t`, `Inf` where nothing was measured.
"""
noise(t::VisTable) = IntensityMap(_noise.(baseimage(weight(t))), domain(t))

_noise(w::Number) = w > 0 ? inv(sqrt(w)) : oftype(w, Inf)
_noise(w::StaticArraysCore.StaticMatrix) = map(_noise, w)

"""
    domain(t::VisTable)

The `(Node, Fr)` domain of `t`. Its coordinates are the baseline `U`, `V` in wavelengths at
each record and channel, and each record's `Ti` (hours since `info(t).epoch`), time-stamp
index `ts` and station indices `ant1`, `ant2`.
"""
domain(t::VisTable) = axisdims(measurement(t))

"""
    info(t::VisTable)

The [`ArrayInfo`](@ref) of the observation that `t` comes from.
"""
info(t::VisTable) = ComradeBase.header(domain(t))

"""
    stations(t::VisTable)

The station grid of `t`: a `DimStack` over `(Ant, Ts)`, the stations and the distinct time
stamps, with the layers `elevation`, `parallactic` and `azimuth` (radians).
"""
stations(t::VisTable) = getfield(t, :stations)

"""
    baselines(t::VisTable)

The station names `(ant1, ant2)` of every record of `t`.
"""
function baselines(t::VisTable)
    names = info(t).stations.name
    d = domain(t)
    return map((i, j) -> (names[i], names[j]), d.ant1, d.ant2)
end

"""
    sites(t::VisTable)

The names of the stations that appear in the records of `t`, in the order of `info(t)`.
"""
function sites(t::VisTable)
    d = domain(t)
    return info(t).stations.name[sort!(unique(vcat(d.ant1, d.ant2)))]
end

Base.size(t::VisTable) = size(measurement(t))
Base.size(t::VisTable, d) = size(measurement(t), d)
Base.length(t::VisTable) = length(measurement(t))

function Base.show(io::IO, ::MIME"text/plain", t::VisTable)
    a = info(t)
    n, nch = size(t)
    kind = eltype(measurement(t)) <: Real ? "amplitudes" : "visibilities"
    println(io, "VisTable: ", kind, ", ", _polname(t.pol))
    println(io, "  source:   ", a.source.name, " (", a.epoch, ")")
    println(io, "  records:  ", n, " × ", nch, " channel", nch == 1 ? "" : "s")
    println(io, "  stations: ", sites(t))
    return print(io, "  measured: ", count(>(0), parent(baseimage(weight(t)))), " of ", length(parent(baseimage(weight(t)))), " values")
end

"""
    getindex(t::VisTable; Node, Fr)

The records and channels of `t` picked by DimensionalData selectors, e.g. `t[Fr = At(230e9)]`
or `t[Node = 1:100]`. A selector that picks a single record or channel keeps its dim. The
station grid is kept whole.
"""
function Base.getindex(t::VisTable; kw...)
    m = measurement(t)
    I = DD.dims2indices(DD.dims(m), DD.Dimensions.kw2dims(values(kw)))
    return _subtable(t, map(i -> i isa Integer ? (i:i) : i, I))
end

_subtable(t::VisTable, I) = VisTable(t.pol, measurement(t)[I...], weight(t)[I...], stations(t))

# The measured values `(record, channel)` of `t`, by record then channel.
_measuredvalues(t::VisTable) = sort!(findall(_measured(t)); by = I -> (I[1], I[2]))

# The time, frequency and station names of each measured value of `t`, in the order of
# `_measuredvalues(t)`.
function _measuredpoints(t::VisTable)
    d = domain(t)
    a = info(t)
    keep = _measuredvalues(t)
    n = getindex.(keep, 1)
    names = a.stations.name
    return (;
        Ti = d.Ti[n], Fr = a.channels.frequency[getindex.(keep, 2)],
        sites = map((i, j) -> (names[i], names[j]), d.ant1[n], d.ant2[n]),
    )
end

# The weight of each value `(record, channel)` of `t`, summed over its components.
_valueweights(t::VisTable) = sum.(baseimage(weight(t)))

# Each station's total data weight, keyed by name.
function _stationweights(t::VisTable)
    names = info(t).stations.name
    d = domain(t)
    w = _valueweights(t)
    total = zeros(length(names))
    for I in CartesianIndices(w)
        n = I[1]
        total[d.ant1[n]] += w[I]
        total[d.ant2[n]] += w[I]
    end
    return NamedTuple{Tuple(names)}(Tuple(total))
end

const _MJD_ZERO = Date(1858, 11, 17)
_mjd(t::VisTable) = Dates.value(info(t).epoch - _MJD_ZERO)

# The numeric storage of a map: `(n, nch)`, or `(n, nch, 2, 2)` for 2×2 values.
_dense(m::IntensityMap{<:SMatrix}) = parent(baseimage(m))
_dense(m::IntensityMap) = baseimage(m)

# Whether each value (record, channel) has at least one measured component.
_measured(t::VisTable) = _ismeasured.(baseimage(weight(t)))
_ismeasured(w) = any(>(0), w)

_withweights(t::VisTable, w) = VisTable(t.pol, measurement(t), IntensityMap(w, domain(t)), stations(t))

"""
    flag(condition, t::VisTable)

`t` with the weights set to `0` wherever `condition(row)` is `true`. A row is one record at one
channel: a `NamedTuple` with the coordinates of that point of `domain(t)` (`U`, `V` in
wavelengths, `Ti`, `Fr`, `ts`, `ant1`, `ant2`), the station names `baseline`, and its
`measurement` and `noise`. Use [`drop`](@ref) to remove the records left
without measurements.

```julia
flag(row -> hypot(row.U, row.V) < 1e8, t)   # flags the short baselines
```
"""
function flag(condition, t::VisTable)
    names = info(t).stations.name
    m, σ = measurement(t), noise(t)
    pts = domainpoints(domain(t))
    keep = map(CartesianIndices(pts)) do I
        p = pts[I]
        row = merge(p, (; baseline = (names[p.ant1], names[p.ant2]), measurement = m[I], noise = σ[I]))
        !condition(row)
    end
    return _withweights(t, baseimage(weight(t)) .* keep)
end

"""
    drop(t::VisTable)

The records of `t` with at least one measured value.
"""
function drop(t::VisTable)
    keep = findall(vec(any(_measured(t); dims = 2)))
    return _subtable(t, (keep, Colon()))
end

"""
    vcat(tables::VisTable...)

One table with the records of all `tables`, e.g. the observations of one source on several
days or in several bands. The tables must hold visibilities of the same polarization and
source. The epoch is the earliest of their epochs, and the times of the other tables are
shifted to it.

The records of one baseline at one time are merged into one record at their mean `u`, `v`.
They must not both measure a channel, and their `u`, `v` must agree to a relative 10⁻³; the
bands of one observation can have coordinates computed a few seconds apart. Join tables
before averaging them: averages over different integrations differ by more.

The stations, channels, scans and tracks are the unions of those of the tables. A station in
several tables must have the same entry in each. A channel a table does not have gets zero
weight in that table's records, and a channel whose center lies inside a different channel is
an error. Overlapping scans are merged into one, and so are overlapping tracks: the bands of
one night share a track, while different nights keep their own.
"""
function Base.vcat(t::VisTable, ts::VisTable...)
    tabs = (t, ts...)
    all(s -> s.pol == t.pol, tabs) ||
        throw(ArgumentError("tables of $(join(unique(map(s -> _polname(s.pol), tabs)), ", ")) cannot be joined"))
    any(s -> eltype(measurement(s)) <: Real, tabs) &&
        throw(ArgumentError("amplitude tables cannot be joined; join the visibility tables, then take amplitudes"))
    source = info(t).source
    all(s -> info(s).source == source, tabs) ||
        throw(ArgumentError("tables of different sources cannot be joined: $(map(s -> info(s).source, tabs))"))
    s = _stack(tabs, _dense ∘ measurement, _dense ∘ weight)
    m = _mergerecords(s.info, s.V, s.W, s.u, s.v, s.Ti, s.ant1, s.ant2, s.rest)
    elevation, parallactic = _joingeometry(s.info, tabs, s.stidx, s.shift, sort!(unique(m.Ti)))
    return VisTable(
        t.pol, m.V; weight = m.W, m.u, m.v, m.Ti, m.ant1, m.ant2, info = s.info, elevation, parallactic,
    )
end

# The earliest epoch of the tables `tabs` and the hours each table's times shift to reach it.
function _epochshift(tabs)
    epoch = minimum(s -> info(s).epoch, tabs)
    return epoch, map(s -> 24.0 * Dates.value(info(s).epoch - epoch), tabs)
end

# The records of `tabs` one after another, with the values `values(s)` and weights `weights(s)`
# of each table `s`, on the stations, channels, scans and tracks of all of them, the time axis
# of the earliest epoch, and the source of the first.
function _stack(tabs, values, weights)
    infos = map(info, tabs)
    epoch, shift = _epochshift(tabs)
    a = ArrayInfo(;
        stations = _joinstations(infos), channels = _joinchannels(tabs),
        scans = _unionintervals(map(a -> a.scans, infos), shift),
        tracks = _unionintervals(map(a -> a.tracks, infos), shift), source = first(infos).source, epoch,
    )

    stidx = map(b -> map(s -> stationindex(a, s), b.stations.name), infos)
    chidx = map(s -> map(f -> findfirst(==(f), a.channels.frequency), _frequencies(s)), tabs)
    Vs, Ws = map(values, tabs), map(weights, tabs)
    fd = size(first(Vs))[3:end]
    all(V -> size(V)[3:end] == fd, Vs) || throw(DimensionMismatch("the tables' values have different shapes"))
    rest = ntuple(_ -> Colon(), length(fd))
    n, nch = sum(s -> size(s, 1), tabs), length(a.channels)
    V = zeros(promote_type(map(eltype, Vs)...), n, nch, fd...)
    W = zeros(promote_type(map(eltype, Ws)...), n, nch, fd...)
    u, v, Ti = (zeros(n) for _ in 1:3)
    ant1, ant2 = zeros(Int, n), zeros(Int, n)
    offset = 0
    for k in eachindex(tabs)
        s = tabs[k]
        r = offset .+ (1:size(s, 1))
        V[r, chidx[k], rest...] .= Vs[k]
        W[r, chidx[k], rest...] .= Ws[k]
        d = domain(s)
        c = ComradeBase.coords(d)
        u[r] .= c.u
        v[r] .= c.v
        Ti[r] .= d.Ti .+ shift[k]
        ant1[r] .= stidx[k][d.ant1]
        ant2[r] .= stidx[k][d.ant2]
        offset += size(s, 1)
    end
    return (; info = a, V, W, u, v, Ti, ant1, ant2, rest, shift, stidx)
end

# A table of the measured values of all `tabs`, unmerged, on the time axis of the earliest epoch,
# each weighted by its total weight: what an instrument's priors are built on. Its values are
# zero, and its `u`, `v` and geometry mix the tables' sources.
function _measuredtable(tabs)
    s = _stack(tabs, t -> zeros(ComplexF64, size(t)), _valueweights)
    return VisTable(IPol(), s.V; weight = s.W, s.u, s.v, s.Ti, s.ant1, s.ant2, info = s.info)
end

const _JOIN_UV_RTOL = 1.0e-3

function _mergerecords(a::ArrayInfo, V, W, u, v, Ti, ant1, ant2, rest)
    records = collect(zip(Ti, ant1, ant2))
    merged = unique(records)
    node = Dict(r => i for (i, r) in pairs(merged))
    nm = length(merged)
    measured = map(I -> any(>(0), view(W, I[1], I[2], rest...)), CartesianIndices((size(W, 1), size(W, 2))))
    seen = falses(nm, size(W, 2))
    Vm, Wm = zeros(eltype(V), nm, size(V)[2:end]...), zeros(eltype(W), nm, size(W)[2:end]...)
    um, vm, count = zeros(nm), zeros(nm), zeros(Int, nm)
    for r in eachindex(records)
        i = node[records[r]]
        any(view(seen, i, :) .& view(measured, r, :)) &&
            throw(ArgumentError("two records of the baseline $(a.stations.name[[ant1[r], ant2[r]]]) at $(Ti[r]) h measure the same channel"))
        seen[i, :] .|= view(measured, r, :)
        selectdim(Vm, 1, i) .+= selectdim(V, 1, r)
        selectdim(Wm, 1, i) .+= selectdim(W, 1, r)
        um[i] += u[r]
        vm[i] += v[r]
        count[i] += 1
    end
    um ./= count
    vm ./= count
    for r in eachindex(records)
        i = node[records[r]]
        hypot(u[r] - um[i], v[r] - vm[i]) <= _JOIN_UV_RTOL * hypot(um[i], vm[i]) ||
            throw(ArgumentError("the records of the baseline $(a.stations.name[[ant1[r], ant2[r]]]) at $(Ti[r]) h have u, v differing by more than a relative $(_JOIN_UV_RTOL); averaged tables can differ more than this, so join tables before averaging them"))
    end
    return (;
        V = Vm, W = Wm, u = um, v = vm,
        Ti = getindex.(merged, 1), ant1 = getindex.(merged, 2), ant2 = getindex.(merged, 3),
    )
end

function _joinstations(infos)
    cols = map(a -> Tables.columntable(a.stations), infos)
    names = keys(first(cols))
    all(c -> issetequal(keys(c), names), cols) ||
        throw(ArgumentError("the station tables have different columns: $(map(keys, cols))"))
    picks = Tuple{Int, Int}[]
    for k in eachindex(cols), i in eachindex(cols[k].name)
        c = cols[k]
        j = findfirst(p -> cols[p[1]].name[p[2]] == c.name[i], picks)
        if isnothing(j)
            push!(picks, (k, i))
        else
            k0, i0 = picks[j]
            all(col -> isequal(cols[k0][col][i0], c[col][i]), names) ||
                throw(ArgumentError("station $(c.name[i]) differs between the tables: $(infos[k0].stations[i0]) and $(infos[k].stations[i])"))
        end
    end
    return NamedTuple{names}(map(col -> [cols[k][col][i] for (k, i) in picks], names))
end

_frequencies(t::VisTable) = collect(DD.lookup(measurement(t), Fr))

# The channels of the `Fr` axis of each table, with their widths from its `info`.
function _joinchannels(tabs)
    frequency = reduce(vcat, map(_frequencies, tabs))
    width = reduce(vcat, map(tabs) do s
        ch = info(s).channels
        map(f -> ch.width[findfirst(==(f), ch.frequency)], _frequencies(s))
    end)
    order = sortperm(frequency)
    f, w = frequency[order], width[order]
    for k in eachindex(f)[2:end]
        if f[k] == f[k - 1]
            w[k] == w[k - 1] ||
                throw(ArgumentError("the channel at $(f[k]) Hz has the widths $(w[k - 1]) and $(w[k]) Hz in different tables"))
        elseif f[k] - f[k - 1] < max(w[k], w[k - 1]) / 2
            throw(ArgumentError("the channels at $(f[k - 1]) and $(f[k]) Hz overlap"))
        end
    end
    keep = [k == firstindex(f) || f[k] != f[k - 1] for k in eachindex(f)]
    return (; frequency = f[keep], width = w[keep])
end

# The union of the intervals of the tables `ivs`, each shifted by its `shift`, with overlapping
# intervals merged.
function _unionintervals(ivs, shift)
    start = reduce(vcat, [iv.start .+ s for (iv, s) in zip(ivs, shift)])
    stop = reduce(vcat, [iv.stop .+ s for (iv, s) in zip(ivs, shift)])
    order = sortperm(start)
    merged = [(start[first(order)], stop[first(order)])]
    for k in order[2:end]
        a, b = merged[end]
        if start[k] <= b
            merged[end] = (a, max(b, stop[k]))
        else
            push!(merged, (start[k], stop[k]))
        end
    end
    return (; start = first.(merged), stop = last.(merged))
end

# The station geometry at the time stamps `times`: computed from the station positions, except
# where a table carries its own.
function _joingeometry(a::ArrayInfo, tabs, stidx, shift, times)
    elevation, parallactic, _ = _geometry(a, times)
    for k in eachindex(tabs)
        g = stations(tabs[k])
        for (j, τ) in pairs(DD.lookup(g, Ts))
            col = searchsortedfirst(times, τ + shift[k])
            col <= lastindex(times) && times[col] == τ + shift[k] || continue
            elevation[stidx[k], col] .= parent(g.elevation)[:, j]
            parallactic[stidx[k], col] .= parent(g.parallactic)[:, j]
        end
    end
    return elevation, parallactic
end

"""
    add_fractional_noise(t::VisTable, f)

`t` with the noise of every measured value raised to `√(σ² + f²|V|²)`. For coherencies, `|V|`
is the Stokes I amplitude from the measured parallel hands.
"""
function add_fractional_noise(t::VisTable, f)
    return _withweights(t, _fractionalweight.(baseimage(weight(t)), _fractionalscale(t), f))
end

_fractionalweight(w::Number, a, f) = w > 0 ? inv(inv(w) + f^2 * a^2) : zero(w)
_fractionalweight(w::StaticArraysCore.StaticMatrix, a, f) = map(wi -> _fractionalweight(wi, a, f), w)

_fractionalscale(t::VisTable) = abs.(baseimage(measurement(t)))
_fractionalscale(t::VisTable{Correlations}) = _stokesiamplitude.(baseimage(measurement(t)), baseimage(weight(t)))

# |I| of one value from its measured parallel hands.
function _stokesiamplitude(v::StaticArraysCore.StaticMatrix, w::StaticArraysCore.StaticMatrix)
    both = w[1, 1] > 0 && w[2, 2] > 0
    return abs(both ? (v[1, 1] + v[2, 2]) / 2 : w[1, 1] > 0 ? v[1, 1] : v[2, 2])
end

"""
    beamsize(t::VisTable)

The inverse of the longest projected baseline, in wavelengths, among the measured values of
`t`: the approximate resolution in radians.
"""
function beamsize(t::VisTable)
    d = domain(t)
    uvd = hypot.(d.U, d.V)
    return inv(maximum(uvd[_measured(t)]))
end

"""
    reset_mounts!(t::VisTable, overrides)
    reset_mounts!(t::VisTable, file::AbstractString)

Updates the station table of `t` in place from a per-station `overrides` dict, as made by
[`load_array_txt`](@ref). Each entry is a `Symbol => NamedTuple` with any of the keys
`fr_parallactic` and `fr_elevation` (the mount becomes `Mount(fr_parallactic, fr_elevation)`)
and `fr_offset` (every feed's receptor angle); other keys are not read.
"""
function reset_mounts!(t::VisTable, overrides)
    st = info(t).stations
    for i in eachindex(st.name)
        ov = get(overrides, st.name[i], nothing)
        isnothing(ov) && continue
        if haskey(ov, :fr_parallactic) || haskey(ov, :fr_elevation)
            m = st.mount[i]
            m isa Mount || throw(ArgumentError("station $(st.name[i]) has a $(nameof(typeof(m))), not a coefficient mount"))
            st.mount[i] = Mount(get(ov, :fr_parallactic, m.parallactic), get(ov, :fr_elevation, m.elevation), m.azimuth, m.offset)
        end
        haskey(ov, :fr_offset) && (st.feed_angles[i] = map(_ -> Float64(ov[:fr_offset]), st.feed_angles[i]))
    end
    return t
end
reset_mounts!(t::VisTable, file::AbstractString) = reset_mounts!(t, load_array_txt(file))

"""
    amplitudes(t::VisTable; debias = false)

The amplitudes `|V|` of the scalar visibilities of `t`, with the same domain and weights. With
`debias = true` they are `√max(|V|² − σ², 0)`, the noise bias removed.
"""
function amplitudes(t::VisTable; debias = false)
    eltype(measurement(t)) <: Complex ||
        throw(ArgumentError("amplitudes are taken of a table of scalar complex visibilities"))
    A = abs.(baseimage(measurement(t)))
    if debias
        w = baseimage(weight(t))
        A = @. ifelse(w > 0, sqrt(max(A^2 - inv(w), zero(A))), A)
    end
    return VisTable(t.pol, IntensityMap(A, domain(t)), weight(t), stations(t))
end
