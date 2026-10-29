module ComradeVLBIFilesExt

using Comrade
using VLBIFiles
using VLBIFiles: VLBIData, VLBI
using Dates
using StaticArraysCore: SMatrix, SVector
using StructArrays: StructArray


# Extract a numeric Hz value from a Unitful Quantity without needing a direct
# Unitful import. `Quantity` has a `.val` field whose units are determined by
# the Quantity's type; for a Hz-typed Quantity this gives Hz.
_to_hz(q) = Float64(q.val)

function _radec(uvd::VLBIFiles.UVData)
    fh = uvd.header.fits
    ra_deg = Float64(VLBIFiles.axis_dict(fh, "RA")["CRVAL"])
    dec_deg = Float64(VLBIFiles.axis_dict(fh, "DEC")["CRVAL"])
    return (ra_deg, dec_deg)
end


# Each time stamp's integration time in seconds.
_inttimes(uvd::VLBIFiles.UVData) = _inttimes(VLBIFiles.uvtable_wide(uvd))
function _inttimes(rows)
    rows = convert(StructArray, rows)
    hasproperty(rows, :int_time) ||
        throw(ArgumentError("the scan boundaries need integration times, and the data give none (no INTTIM column)"))
    Δ = Dict{DateTime, Float64}()
    for (t, d) in zip(rows.datetime, rows.int_time)
        Δ[t] = max(get(Δ, t, 0.0), Float64(d.val))
    end
    return Δ
end

# The scans of `rows`, by their scan ids or else `GapBasedScans`, each from the start of its
# first integration to the end of its last.
function _scan_intervals(rows, inttime = _inttimes(rows))
    rows = convert(StructArray, rows)
    ids = hasproperty(rows, :scan_id) ? rows.scan_id : VLBIData.scan_ids(VLBI.GapBasedScans(), rows)
    starts, stops = Dict{eltype(ids), DateTime}(), Dict{eltype(ids), DateTime}()
    for (k, t) in zip(ids, rows.datetime)
        half = Millisecond(round(Int, 500 * inttime[t]))
        starts[k] = min(get(starts, k, t - half), t - half)
        stops[k] = max(get(stops, k, t + half), t + half)
    end
    return [VLBIData.IntervalSets.Interval(starts[k], stops[k]) for k in sort!(collect(keys(starts)))]
end

# VisTable -------------------------------------------------------------------

_hours(dt, epoch) = Dates.value(dt - DateTime(epoch)) / 3.6e6

const _MOUNTS = Dict(
    VLBIData.AntennaMountType.AltAzimuth => MountAltAz(),
    VLBIData.AntennaMountType.Equatorial => MountEquatorial(),
    VLBIData.AntennaMountType.Orbiting => MountOrbiting(),
    VLBIData.AntennaMountType.XY => MountXY(),
    VLBIData.AntennaMountType.NaismithR => MountNasmythR(),
    VLBIData.AntennaMountType.NaismithL => MountNasmythL(),
)
# An aperture array or a mount the file does not name has no known feed rotation.
_mount(ant) = get(_MOUNTS, ant.mount_type, MountOther())

const _FEEDS = Dict(:R => RPol(), :L => LPol(), :X => XPol(), :Y => YPol())
function _feeds(ant)
    all(p -> haskey(_FEEDS, p), ant.poltypes) ||
        throw(ArgumentError("station $(ant.name) has the feed polarizations $(ant.poltypes); expected R, L, X or Y"))
    return map(p -> _FEEDS[p], ant.poltypes)
end

function _arrayinfo(antarray, source, channels, epoch, sint)
    ants = unique(a -> a.name, [ant for (_, ant) in pairs(antarray.antennas)])
    stations = (;
        name = [a.name for a in ants],
        position = [SVector{3, Float64}(a.xyz) for a in ants],
        mount = [_mount(a) for a in ants],
        feeds = [_feeds(a) for a in ants],
        # POLAA for both feeds: writers differ on whether POLAB is the second feed's own offset
        # or its dipole angle, POLAA + 90°. NaN where the file gives none; a polarized fit then
        # stops at that station.
        feed_angles = [ntuple(_ -> Float64(first(a.feed_offsets)), 2) for a in ants],
    )
    ends = VLBIData.IntervalSets.endpoints.(sint)
    scans = (; start = [_hours(first(e), epoch) for e in ends], stop = [_hours(last(e), epoch) for e in ends])
    return Comrade.ArrayInfo(;
        stations, channels, scans, source = Comrade.Source(source.name, source.ra, source.dec), epoch,
    )
end

# Every channel's center frequency and width in Hz, in the file's order.
function _channels(fws)
    frequency = map(_to_hz, collect(VLBIFiles.frequencies(fws)))
    width = reduce(vcat, [fill(abs(_to_hz(fw.width)) / fw.nchan, fw.nchan) for fw in fws])
    return (; frequency, width)
end

# Where `pol` sits among the correlation products `products` of a record: the parallel hands
# for Stokes I, the feed-product slots for a matrix.
function _productlayout(::IPol, products)
    hands = findall(p -> p === :I || p in (:RR, :LL, :XX, :YY), products)
    length(hands) in (1, 2) || throw(ArgumentError("expected one or two parallel-hand products, got $products"))
    return (first(hands), last(hands))
end
_productlayout(pol, products) = VLBIData.feed_slot_map(products)
_layouttype(::IPol) = NTuple{2, Int}
_layouttype(pol) = SMatrix{2, 2, Int8, 4}

# The value and noise of `pol` from a record's visibilities `vis` and weights `w` at channel `k`,
# with `NaN` and `Inf` where nothing is measured. Stokes I is the mean of the measured parallel
# hands.
function _recordvalue(::IPol, hands, vis, w, k)
    n, total, var = 0, zero(ComplexF64), 0.0
    for (j, h) in pairs(hands)
        j == 2 && h == hands[1] && continue
        wh = Float64(w[h, k])
        wh > 0 || continue
        n += 1
        total += ComplexF64(vis[h, k])
        var += inv(wh)
    end
    return iszero(n) ? (complex(NaN), Inf) : (total / n, sqrt(var) / n)
end
function _recordvalue(pol, slots, vis, w, k)
    m = map(i -> i > 0 && w[i, k] > 0 ? ComplexF64(vis[i, k]) : complex(NaN), slots)
    s = map(i -> i > 0 && w[i, k] > 0 ? inv(sqrt(Float64(w[i, k]))) : Inf, slots)
    return _select(pol, SMatrix{2, 2, ComplexF64, 4}(m), SMatrix{2, 2, Float64, 4}(s))
end
_select(::Comrade.Correlations, m, s) = (m, s)
_select(::Comrade.Correlation{I, J}, m, s) where {I, J} = (m[I, J], s[I, J])

"""
    VisTable(uvd::VLBIFiles.UVData, pol)

Every integration and channel of a uvfits file loaded with VLBIFiles, as a `VisTable` of the
polarization `pol`: `IPol()`, a correlation such as `Correlation(1, 1)`, or `Correlations()`
(see [`VisTable`](@ref Comrade.VisTable)). Stokes I is the mean of the measured parallel hands.
Average the result with [`average`](@ref Comrade.average), e.g. `average(t, ByScan(), ChannelBins(4))`.
"""
function Comrade.VisTable(uvd::VLBIFiles.UVData, pol::Comrade.TablePolarization)
    wide = VLBIFiles.uvtable_wide(uvd)
    chans = _channels(uvd.freq_windows)
    order = sortperm(chans.frequency)
    channels = map(c -> c[order], chans)
    allunique(channels.frequency) || throw(ArgumentError("the file has two channels at one frequency"))
    column = invperm(order)

    fd = Comrade._valueshape(pol)
    n, nch = length(wide), length(order)
    V = fill(complex(NaN), n, nch, fd...)
    σ = fill(Inf, n, nch, fd...)
    # The file's product labels hold for every baseline when all stations have the same feeds.
    stokes = uvd.header.stokes
    asis = length(unique(a.poltypes for arr in uvd.ant_arrays for a in arr.antennas)) == 1
    layouts = Dict{eltype(wide.baseline), _layouttype(pol)}()
    for i in eachindex(wide.baseline, wide.visibility, wide.weight)
        bl = wide.baseline[i]
        layout = get!(layouts, bl) do
            _productlayout(pol, asis ? stokes : map(s -> VLBIFiles.to_proper_stokes(s, uvd.ant_arrays, bl), stokes))
        end
        vis, w = wide.visibility[i], wide.weight[i]
        for k in eachindex(column)
            m, s = _recordvalue(pol, layout, vis, w, k)
            selectdim(selectdim(V, 1, i), 1, column[k]) .= m
            selectdim(selectdim(σ, 1, i), 1, column[k]) .= s
        end
    end
    keep = findall(i -> any(isfinite, selectdim(σ, 1, i)), axes(σ, 1))
    isempty(keep) && throw(ArgumentError("the file has no measured values of $pol"))
    recordkeys = [(wide.datetime[i], Tuple(wide.baseline[i].antennas)) for i in keep]
    allunique(recordkeys) || throw(ArgumentError("the file has two records of one baseline at one time"))

    epoch = Date(minimum(wide.datetime[keep]))
    sint = _scan_intervals(StructArray((; datetime = wide.datetime[keep])), _inttimes(wide))
    ra, dec = _radec(uvd)
    info = _arrayinfo(only(uvd.ant_arrays), (; name = uvd.header.object, ra, dec), channels, epoch, sint)
    uvw = wide.uvw[keep]
    return Comrade.VisTable(
        pol, V[keep, ntuple(_ -> Colon(), ndims(V) - 1)...], σ[keep, ntuple(_ -> Colon(), ndims(σ) - 1)...];
        u = [Float64(p[1].val) for p in uvw], v = [Float64(p[2].val) for p in uvw],
        Ti = [_hours(k[1], epoch) for k in recordkeys],
        ant1 = [k[2][1] for k in recordkeys], ant2 = [k[2][2] for k in recordkeys], info,
    )
end

end # module
