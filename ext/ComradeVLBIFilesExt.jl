module ComradeVLBIFilesExt

using Comrade
using VLBIFiles
using VLBIData
using VLBIData: VLBI, AntennaMountType, uvtable
using Comrade: PolarizedTypes
using Comrade.PolarizedTypes: CirBasis, LinBasis
using AstroLib
using Dates
using LinearAlgebra
using StaticArraysCore: SMatrix
using StructArrays: StructVector, StructArray


# Extract a numeric Hz value from a Unitful Quantity without needing a direct
# Unitful import. `Quantity` has a `.val` field whose units are determined by
# the Quantity's type; for a Hz-typed Quantity this gives Hz.
_to_hz(q) = Float64(q.val)


# Mount-type → (fr_parallactic, fr_elevation, fr_offset_rad).
# Conventions match eht-imaging's antenna table: 1 = field rotates with the
# corresponding angle, 0 = no rotation, ±1 for Naismith handedness.
const _MOUNT_FR = Dict(
    AntennaMountType.AltAzimuth => (1.0, 0.0, 0.0),
    AntennaMountType.Equatorial => (0.0, 0.0, 0.0),
    AntennaMountType.Orbiting => (0.0, 0.0, 0.0),
    AntennaMountType.XY => (1.0, 0.0, 0.0),
    AntennaMountType.NaismithR => (1.0, 1.0, 0.0),
    AntennaMountType.NaismithL => (1.0, -1.0, 0.0),
    AntennaMountType.ApertureArray => (0.0, 0.0, 0.0),
    AntennaMountType.Unknown => (1.0, 0.0, 0.0),
)


function _mount_coeffs(ant, mount_overrides)
    if mount_overrides !== nothing && haskey(mount_overrides, ant.name)
        return mount_overrides[ant.name]
    end
    return get(_MOUNT_FR, ant.mount_type, (1.0, 0.0, 0.0))
end


function _hour_of_day(dt::DateTime)
    midnight = DateTime(Date(dt))
    ms = Dates.value(dt - midnight)        # milliseconds since midnight
    return ms / 3.6e6
end


function _jd(dt::DateTime)
    # Julian date: JD 2440587.5 = 1970-01-01T00:00:00 UTC
    s = Dates.datetime2unix(dt)
    return 2440587.5 + s / 86400
end


function _radec(uvd::VLBIFiles.UVData)
    fh = uvd.header.fits
    ra_deg = Float64(VLBIFiles.axis_dict(fh, "RA")["CRVAL"])
    dec_deg = Float64(VLBIFiles.axis_dict(fh, "DEC")["CRVAL"])
    return (ra_deg, dec_deg)
end


function _mjd(uvd::VLBIFiles.UVData)
    return Dates.value(uvd.header.date_obs - Date(1858, 11, 17))
end


function _build_tarr(
        uvd::VLBIFiles.UVData; mount_overrides = nothing,
        array_overrides = nothing
    )
    # collect unique antennas across (typically just one) ant_array
    seen = Set{Symbol}()
    sites = Symbol[]
    X = Float64[]; Y = Float64[]; Z = Float64[]
    SEFD1 = Float64[]; SEFD2 = Float64[]
    fr_par = Float64[]; fr_el = Float64[]; fr_off = Float64[]
    for ar in uvd.ant_arrays
        for (_, ant) in pairs(ar.antennas)
            ant.name in seen && continue
            push!(seen, ant.name)
            push!(sites, ant.name)
            push!(X, ant.xyz[1])
            push!(Y, ant.xyz[2])
            push!(Z, ant.xyz[3])
            push!(SEFD1, NaN)
            push!(SEFD2, NaN)
            fp, fe, fo = _mount_coeffs(ant, mount_overrides)
            push!(fr_par, fp); push!(fr_el, fe); push!(fr_off, fo)
        end
    end
    if array_overrides !== nothing
        for i in eachindex(sites)
            ov = get(array_overrides, sites[i], nothing)
            ov === nothing && continue
            haskey(ov, :SEFD1)          && (SEFD1[i] = ov[:SEFD1])
            haskey(ov, :SEFD2)          && (SEFD2[i] = ov[:SEFD2])
            haskey(ov, :fr_parallactic) && (fr_par[i] = ov[:fr_parallactic])
            haskey(ov, :fr_elevation)   && (fr_el[i] = ov[:fr_elevation])
            haskey(ov, :fr_offset)      && (fr_off[i] = ov[:fr_offset])
        end
    end
    return StructArray(
        (;
            sites, X, Y, Z, SEFD1, SEFD2,
            fr_parallactic = fr_par,
            fr_elevation = fr_el,
            fr_offset = fr_off,
        )
    )
end


"""
    load_array_txt(path) -> Dict{Symbol, NamedTuple}

Read an ehtim-style antenna text file and return a dictionary keyed by site
symbol. Each value is a NamedTuple with the columns present in the file:
`(SEFD1, SEFD2, fr_parallactic, fr_elevation, fr_offset)` plus the d-term
fields `(DR, DL)` if present. `fr_offset` is converted from degrees (the
ehtim file convention) to radians.

The parsed file is in the same column layout that `eht-imaging` uses:

    SITE  X  Y  Z  SEFDR  SEFDL  DR_re DR_im  DL_re DL_im  FR_PAR  FR_ELEV  FR_OFF

`X, Y, Z` are ignored here — antenna positions come from the uvfits AIPS-AN
table. The result can be passed to any `extract_*` as
`array_overrides=load_array_txt("array.txt")`.
"""
function Comrade.load_array_txt(path::AbstractString)
    overrides = Dict{Symbol, NamedTuple}()
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        toks = split(line)
        length(toks) >= 13 || continue
        site = Symbol(toks[1])
        sefd_r = parse(Float64, toks[5])
        sefd_l = parse(Float64, toks[6])
        dr = complex(parse(Float64, toks[7]), parse(Float64, toks[8]))
        dl = complex(parse(Float64, toks[9]), parse(Float64, toks[10]))
        fr_par = parse(Float64, toks[11])
        fr_el = parse(Float64, toks[12])
        fr_off = deg2rad(parse(Float64, toks[13]))
        overrides[site] = (
            SEFD1 = sefd_r, SEFD2 = sefd_l,
            fr_parallactic = fr_par, fr_elevation = fr_el, fr_offset = fr_off,
            DR = dr, DL = dl,
        )
    end
    return overrides
end


function _build_scans(uvtbl)
    times = sort(unique(uvtbl.datetime))
    isempty(times) && return StructArray((; start = Float64[], stop = Float64[]))
    day0 = Date(times[1])
    hrs = map(t -> _hour_of_day(t) + 24 * Dates.value(Date(t) - day0), times)
    return StructArray((; start = hrs, stop = hrs))
end


function _arrayconfig(
        uvd::VLBIFiles.UVData, uvtbl;
        mount_overrides = nothing, array_overrides = nothing
    )
    ra_deg, dec_deg = _radec(uvd)
    # Pyehtim/ehtim stores RA in hours; match that convention so downstream
    # tooling that consumes either path agrees on units.
    ra_hours = ra_deg / 15
    mjd = _mjd(uvd)
    source = Symbol(uvd.header.object)
    bw = sum(fw -> _to_hz(fw.width), uvd.freq_windows)
    tarr = _build_tarr(uvd; mount_overrides, array_overrides)
    scans = _build_scans(uvtbl)

    n = length(uvtbl)
    U = Vector{Float64}(undef, n)
    V = Vector{Float64}(undef, n)
    Ti = Vector{Float64}(undef, n)
    Fr = Vector{Float64}(undef, n)
    sites_v = Vector{Tuple{Symbol, Symbol}}(undef, n)
    elevation = Vector{Tuple{Float64, Float64}}(undef, n)
    parallactic = Vector{Tuple{Float64, Float64}}(undef, n)

    # caches built from antenna metadata
    site_xyz = Dict{Symbol, NTuple{3, Float64}}()
    poltype_map = Dict{Symbol, NTuple{2, Symbol}}()
    for ar in uvd.ant_arrays, (_, ant) in pairs(ar.antennas)
        site_xyz[ant.name] = (Float64(ant.xyz[1]), Float64(ant.xyz[2]), Float64(ant.xyz[3]))
        poltype_map[ant.name] = ant.poltypes
    end

    # Build per-row polbasis. We assume all stations in this file share a single
    # polbasis (matching ehtim semantics); error out if they don't.
    pol_set = unique(values(poltype_map))
    length(pol_set) == 1 || error("Mixed polarization bases across stations not yet supported: $pol_set")
    pb_first = pol_set[1][1]
    single_pb = pb_first in (:R, :L) ? (CirBasis(), CirBasis()) : (LinBasis(), LinBasis())
    polbasis = fill(single_pb, n)

    geocache = Dict{Tuple{Symbol, DateTime}, NTuple{2, Float64}}()

    @inbounds for i in 1:n
        r = uvtbl[i]
        s1, s2 = r.spec.bl.antennas
        U[i] = Float64(r.spec.uv.u)
        V[i] = Float64(r.spec.uv.v)
        Ti[i] = _hour_of_day(r.datetime)
        Fr[i] = _to_hz(r.freq_spec.freq)
        sites_v[i] = (s1, s2)

        jd = _jd(r.datetime)
        e1, p1 = get!(geocache, (s1, r.datetime)) do
            Comrade.elevation_parallactic(site_xyz[s1], ra_deg, dec_deg, jd)
        end
        e2, p2 = get!(geocache, (s2, r.datetime)) do
            Comrade.elevation_parallactic(site_xyz[s2], ra_deg, dec_deg, jd)
        end
        elevation[i] = (e1, e2)
        parallactic[i] = (p1, p2)
    end

    data = StructArray{Comrade.EHTArrayBaselineDatum{Float64, eltype(polbasis), Float64}}(
        (; U, V, Ti, Fr, sites = sites_v, polbasis, elevation, parallactic)
    )
    return Comrade.EHTArrayConfiguration(
        bw, tarr, scans, mjd, ra_hours, dec_deg,
        source, :UTC, data
    )
end


# polarization-row filtering helpers -------------------------------------

# Map Comrade pol Symbol -> the set of uvfits stokes labels and the linear
# combination weights to apply to (RR, LL, RL, LR) or (XX, YY, XY, YX) rows.
# Returns nothing if `pol` is one of the raw stokes labels (use directly).
function _pol_filter(pol::Symbol, available::Set{Symbol})
    if pol in available
        return :raw
    end
    pol === :I && (:RR in available && :LL in available) && return :I_circ
    pol === :I && (:XX in available && :YY in available) && return :I_lin
    error("pol=$pol not available; uvtable stokes are $available")
end


# Load the uvtable, optionally apply averaging, then filter to the requested
# polarization. The averaging strategies come straight from VLBIData:
#   frequency_average=true (or anything truthy) -> VLBI.average_data(VLBI.ByFrequency(), ...)
#   time_average=VLBI.GapBasedScans()           -> VLBI.average_data(time_average, ...)
#   time_average=VLBI.FixedTimeIntervals(...)
function _prep_uvtable(
        uvd::VLBIFiles.UVData; pol::Symbol,
        time_average = nothing, frequency_average = true
    )
    rows = VLBIData.uvtable(uvd)
    if frequency_average !== nothing && frequency_average !== false
        rows = VLBI.average_data(VLBI.ByFrequency(), rows)
    end
    if time_average !== nothing && time_average !== false
        rows = VLBI.average_data(time_average, rows)
    end
    return _filter_to_pol(rows, pol)
end


function _filter_to_pol(uvtbl, pol::Symbol)
    available = Set(unique(uvtbl.stokes))
    mode = _pol_filter(pol, available)
    if mode === :raw
        return uvtbl[uvtbl.stokes .== pol]
    elseif mode === :I_circ
        # take I = (RR + LL) / 2 per (datetime, baseline, freq_spec)
        return _combine_parallel_hands(uvtbl, (:RR, :LL))
    elseif mode === :I_lin
        return _combine_parallel_hands(uvtbl, (:XX, :YY))
    end
end


# Average two parallel-hand stokes per group into a single Stokes-I row.
# A row only appears in the output if BOTH hands are present for that group.
function _combine_parallel_hands(uvtbl, hands::NTuple{2, Symbol})
    h1, h2 = hands
    rows1 = uvtbl[uvtbl.stokes .== h1]
    rows2 = uvtbl[uvtbl.stokes .== h2]
    keyfn = r -> (r.datetime, r.spec.bl.antennas, r.freq_spec.freq, r.spec.uv.u, r.spec.uv.v)
    map2 = Dict{Any, Int}()
    for (j, r) in pairs(rows2)
        map2[keyfn(r)] = j
    end

    keep1 = Int[]
    keep2 = Int[]
    for (i, r) in pairs(rows1)
        j = get(map2, keyfn(r), 0)
        j == 0 && continue
        push!(keep1, i)
        push!(keep2, j)
    end
    r1 = rows1[keep1]
    r2 = rows2[keep2]

    # construct Stokes-I rows: I = (h1 + h2)/2; sigma_I = 0.5 * sqrt(σ1² + σ2²)
    n = length(r1)
    new_value = Vector{eltype(r1.value)}(undef, n)
    for i in 1:n
        v1, e1 = r1[i].value.v, r1[i].value.u
        v2, e2 = r2[i].value.v, r2[i].value.u
        I = (v1 + v2) / 2
        σI = sqrt(e1^2 + e2^2) / 2
        new_value[i] = typeof(r1[i].value)(I, σI)
    end

    return StructArray(
        (;
            datetime = r1.datetime,
            stokes = fill(:I, n),
            freq_spec = r1.freq_spec,
            spec = r1.spec,
            value = new_value,
        )
    )
end


function getvisfield(uvtbl)
    n = length(uvtbl)
    vis = Vector{ComplexF64}(undef, n)
    err = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        vis[i] = ComplexF64(uvtbl[i].value.v)
        err[i] = Float64(uvtbl[i].value.u)
    end
    return vis, err
end


# extract_vis ---------------------------------------------------------------

"""
    extract_vis(uvd::VLBIFiles.UVData; pol=:I, mount_overrides=nothing, kwargs...)

Build a Comrade `EHTObservationTable` of complex visibilities from a VLBIFiles
`UVData` object. No averaging is performed: the output has one row for each
row of `VLBIFiles.uvtable(uvd)` after polarization filtering.

`pol`:
- `:I` (default) — combine parallel hands (RR+LL or XX+YY) into Stokes I.
- A raw stokes label like `:RR`, `:LL`, `:XX`, `:YY` — keep only matching rows.

`mount_overrides::Dict{Symbol, NTuple{3, Float64}}` lets you override the
`(fr_parallactic, fr_elevation, fr_offset)` triplet derived from the uvfits
mount type for individual stations.
"""
function Comrade.extract_vis(
        uvd::VLBIFiles.UVData; pol = :I,
        time_average = nothing,
        frequency_average = true,
        mount_overrides = nothing,
        array_overrides = nothing,
        arrayfile = nothing, kwargs...
    )
    array_overrides === nothing && arrayfile !== nothing &&
        (array_overrides = Comrade.load_array_txt(arrayfile))
    rows = _prep_uvtable(uvd; pol, time_average, frequency_average)
    config = _arrayconfig(uvd, rows; mount_overrides, array_overrides)
    vis, err = getvisfield(rows)
    T = Comrade.EHTVisibilityDatum{pol, eltype(err), typeof(config[1])}
    return Comrade.EHTObservationTable{T}(vis, err, config)
end


"""
    extract_amp(uvd::VLBIFiles.UVData; pol=:I, debias=false,
                 mount_overrides=nothing, kwargs...)

Build a Comrade `EHTObservationTable` of visibility amplitudes. If
`debias=true` the Rice-bias correction `√max(|V|² − σ², 0)` is applied.
"""
function Comrade.extract_amp(
        uvd::VLBIFiles.UVData; pol = :I, debias = false,
        time_average = nothing,
        frequency_average = true,
        mount_overrides = nothing,
        array_overrides = nothing,
        arrayfile = nothing, kwargs...
    )
    array_overrides === nothing && arrayfile !== nothing &&
        (array_overrides = Comrade.load_array_txt(arrayfile))
    rows = _prep_uvtable(uvd; pol, time_average, frequency_average)
    config = _arrayconfig(uvd, rows; mount_overrides, array_overrides)
    vis, err = getvisfield(rows)
    amp = abs.(vis)
    if debias
        amp = @. sqrt(max(amp^2 - err^2, zero(amp)))
    end
    T = Comrade.EHTVisibilityAmplitudeDatum{pol, eltype(amp), typeof(config[1])}
    return Comrade.EHTObservationTable{T}(amp, err, config)
end


# extract_coherency ---------------------------------------------------------

"""
    extract_coherency(uvd::VLBIFiles.UVData; mount_overrides=nothing, kwargs...)

Build a Comrade `EHTObservationTable` of 2×2 coherency matrices. Requires the
uvfits to contain all four circular hands (RR, LL, RL, LR). Missing hands for
a given (time, baseline, frequency) row are filled with `NaN+NaN*im` and
infinite uncertainty.
"""
function Comrade.extract_coherency(
        uvd::VLBIFiles.UVData;
        time_average = nothing,
        frequency_average = true,
        mount_overrides = nothing,
        array_overrides = nothing,
        arrayfile = nothing, kwargs...
    )
    array_overrides === nothing && arrayfile !== nothing &&
        (array_overrides = Comrade.load_array_txt(arrayfile))
    uvtbl_full = VLBIData.uvtable(uvd)
    if frequency_average !== nothing
        uvtbl_full = VLBI.average_data(VLBI.ByFrequency(), uvtbl_full)
    end
    if time_average !== nothing
        uvtbl_full = VLBI.average_data(time_average, uvtbl_full)
    end
    available = Set(unique(uvtbl_full.stokes))
    needed = (:RR, :RL, :LR, :LL)
    for h in needed
        h in available || error("uvfits is missing $h hand; coherency requires all four circular hands")
    end

    keyfn = r -> (r.datetime, r.spec.bl.antennas, r.freq_spec.freq, r.spec.uv.u, r.spec.uv.v)
    rows_by = Dict{Symbol, Any}()
    for h in needed
        rows_by[h] = uvtbl_full[uvtbl_full.stokes .== h]
    end

    keys_for = Dict{Symbol, Dict{Any, Int}}()
    for h in needed
        d = Dict{Any, Int}()
        rh = rows_by[h]
        for (i, r) in pairs(rh)
            d[keyfn(r)] = i
        end
        keys_for[h] = d
    end

    # use RR rows as the row "spine"; ensure all four hands have entries
    spine = rows_by[:RR]
    n = length(spine)
    rrv = Vector{ComplexF64}(undef, n); rre = Vector{Float64}(undef, n)
    rlv = Vector{ComplexF64}(undef, n); rle = Vector{Float64}(undef, n)
    lrv = Vector{ComplexF64}(undef, n); lre = Vector{Float64}(undef, n)
    llv = Vector{ComplexF64}(undef, n); lle = Vector{Float64}(undef, n)

    keep = Int[]
    for (i, r) in pairs(spine)
        k = keyfn(r)
        jrl = get(keys_for[:RL], k, 0)
        jlr = get(keys_for[:LR], k, 0)
        jll = get(keys_for[:LL], k, 0)
        # require LL present at minimum; missing crosshands -> NaN
        jll == 0 && continue
        push!(keep, i)
        idx = length(keep)
        rrv[idx] = ComplexF64(r.value.v); rre[idx] = Float64(r.value.u)
        if jrl == 0
            rlv[idx] = ComplexF64(NaN, NaN); rle[idx] = Inf
        else
            rrl = rows_by[:RL][jrl]
            rlv[idx] = ComplexF64(rrl.value.v); rle[idx] = Float64(rrl.value.u)
        end
        if jlr == 0
            lrv[idx] = ComplexF64(NaN, NaN); lre[idx] = Inf
        else
            rlr = rows_by[:LR][jlr]
            lrv[idx] = ComplexF64(rlr.value.v); lre[idx] = Float64(rlr.value.u)
        end
        rll = rows_by[:LL][jll]
        llv[idx] = ComplexF64(rll.value.v); lle[idx] = Float64(rll.value.u)
    end
    nkeep = length(keep)
    resize!(rrv, nkeep); resize!(rre, nkeep)
    resize!(rlv, nkeep); resize!(rle, nkeep)
    resize!(lrv, nkeep); resize!(lre, nkeep)
    resize!(llv, nkeep); resize!(lle, nkeep)

    spine_kept = spine[keep]
    config = _arrayconfig(uvd, spine_kept; mount_overrides, array_overrides)

    cohmat = StructArray{SMatrix{2, 2, ComplexF64, 4}}((rrv, lrv, rlv, llv))
    errmat = StructArray{SMatrix{2, 2, Float64, 4}}((rre, lre, rle, lle))
    T = Comrade.EHTCoherencyDatum{Float64, typeof(config[1]), eltype(cohmat), eltype(errmat)}
    return Comrade.EHTObservationTable{T}(cohmat, errmat, config)
end


# extract_cphase / extract_lcamp -------------------------------------------

# Build the Comrade closure-spec StructArray from a visibility EHTObservationTable
# by enumerating maximal triangles/quadrangles within each scan.
function _build_closures(dvis::Comrade.EHTObservationTable, ::Val{:cphase})
    st = Comrade.timetable(dvis)
    Tvec = Float64[]
    Fvec = Float64[]
    nvec = Float64[]
    bls = Tuple{Symbol, Symbol, Symbol}[]
    for si in eachindex(st)
        s = st[si]
        scan = s.scan
        ant1, ant2 = Comrade.baseline(s)
        meas = Comrade.measurement(scan)
        nz = Comrade.noise(scan)
        # build baseline -> (idx) map for this scan
        blmap = Dict{Tuple{Symbol, Symbol}, Int}()
        for k in eachindex(ant1)
            blmap[(ant1[k], ant2[k])] = k
            blmap[(ant2[k], ant1[k])] = -k
        end
        scan_sites = sort(unique(vcat(ant1, ant2)))
        nsites = length(scan_sites)
        nsites < 3 && continue
        for ii in 1:(nsites - 2), jj in (ii + 1):(nsites - 1), kk in (jj + 1):nsites
            a, b, c = scan_sites[ii], scan_sites[jj], scan_sites[kk]
            i_ab = get(blmap, (a, b), 0)
            i_bc = get(blmap, (b, c), 0)
            i_ca = get(blmap, (c, a), 0)
            (i_ab == 0 || i_bc == 0 || i_ca == 0) && continue
            v_ab = i_ab > 0 ? meas[i_ab] : conj(meas[-i_ab])
            v_bc = i_bc > 0 ? meas[i_bc] : conj(meas[-i_bc])
            v_ca = i_ca > 0 ? meas[i_ca] : conj(meas[-i_ca])
            σ_ab = nz[abs(i_ab)]
            σ_bc = nz[abs(i_bc)]
            σ_ca = nz[abs(i_ca)]
            # closure-phase noise (small-error approx)
            cpnoise = sqrt(
                (σ_ab / abs(v_ab))^2 +
                    (σ_bc / abs(v_bc))^2 +
                    (σ_ca / abs(v_ca))^2
            )
            push!(Tvec, s.time)
            push!(Fvec, scan[1].baseline.Fr)
            push!(nvec, cpnoise)
            push!(bls, (a, b, c))
        end
    end
    return StructArray((; T = Tvec, F = Fvec, noise = nvec, baseline = bls))
end


function _build_closures(dvis::Comrade.EHTObservationTable, ::Val{:lcamp})
    st = Comrade.timetable(dvis)
    Tvec = Float64[]
    Fvec = Float64[]
    nvec = Float64[]
    bls = Tuple{Symbol, Symbol, Symbol, Symbol}[]
    for si in eachindex(st)
        s = st[si]
        scan = s.scan
        ant1, ant2 = Comrade.baseline(s)
        meas = Comrade.measurement(scan)
        nz = Comrade.noise(scan)
        blmap = Dict{Tuple{Symbol, Symbol}, Int}()
        for k in eachindex(ant1)
            blmap[(ant1[k], ant2[k])] = k
            blmap[(ant2[k], ant1[k])] = -k
        end
        scan_sites = sort(unique(vcat(ant1, ant2)))
        nsites = length(scan_sites)
        nsites < 4 && continue
        for ia in 1:(nsites - 3), ib in (ia + 1):(nsites - 2), ic in (ib + 1):(nsites - 1), id in (ic + 1):nsites
            a, b, c, d = scan_sites[ia], scan_sites[ib], scan_sites[ic], scan_sites[id]
            i_ab = get(blmap, (a, b), 0)
            i_cd = get(blmap, (c, d), 0)
            i_ad = get(blmap, (a, d), 0)
            i_bc = get(blmap, (b, c), 0)
            (i_ab == 0 || i_cd == 0 || i_ad == 0 || i_bc == 0) && continue
            v_ab = abs(meas[abs(i_ab)])
            v_cd = abs(meas[abs(i_cd)])
            v_ad = abs(meas[abs(i_ad)])
            v_bc = abs(meas[abs(i_bc)])
            σ_ab = nz[abs(i_ab)]
            σ_cd = nz[abs(i_cd)]
            σ_ad = nz[abs(i_ad)]
            σ_bc = nz[abs(i_bc)]
            lcanoise = sqrt(
                (σ_ab / v_ab)^2 + (σ_cd / v_cd)^2 +
                    (σ_ad / v_ad)^2 + (σ_bc / v_bc)^2
            )
            push!(Tvec, s.time)
            push!(Fvec, scan[1].baseline.Fr)
            push!(nvec, lcanoise)
            push!(bls, (a, b, c, d))
        end
    end
    return StructArray((; T = Tvec, F = Fvec, baseline = bls, noise = nvec))
end


"""
    extract_cphase(uvd::VLBIFiles.UVData; pol=:I, count="min",
                    mount_overrides=nothing, kwargs...)

Extract closure phases. Closure triangles are enumerated within each scan and
then optionally reduced to a minimal independent set (`count="min"`) or kept
as the maximal set (`count="max"`).
"""
function Comrade.extract_cphase(
        uvd::VLBIFiles.UVData; pol = :I, count = "min",
        time_average = nothing,
        frequency_average = true,
        mount_overrides = nothing, kwargs...
    )
    dvis = Comrade.extract_vis(
        uvd; pol, time_average, frequency_average,
        mount_overrides, kwargs...
    )
    cphase = _build_closures(dvis, Val(:cphase))
    clac = Comrade.build_closure_config(dvis, cphase; type = :cphase, count)
    T = Comrade.EHTClosurePhaseDatum{pol, eltype(cphase.T), typeof(arrayconfig(dvis)[1])}
    cp = Comrade.closure_phases(measurement(dvis), Comrade.designmat(clac))
    cp_sig = abs2.(noise(dvis) ./ measurement(dvis))
    cp_cov = Comrade.designmat(clac) * Diagonal(cp_sig) * transpose(Comrade.designmat(clac))
    return Comrade.EHTObservationTable{T}(cp, cp_cov, clac)
end


"""
    extract_lcamp(uvd::VLBIFiles.UVData; pol=:I, count="min",
                   mount_overrides=nothing, kwargs...)

Extract log-closure amplitudes.
"""
function Comrade.extract_lcamp(
        uvd::VLBIFiles.UVData; pol = :I, count = "min",
        time_average = nothing,
        frequency_average = true,
        mount_overrides = nothing, kwargs...
    )
    dvis = Comrade.extract_vis(
        uvd; pol, time_average, frequency_average,
        mount_overrides, kwargs...
    )
    lcamp = _build_closures(dvis, Val(:lcamp))
    clac = Comrade.build_closure_config(dvis, lcamp; type = :lcamp, count)
    cldmat = Comrade.designmat(clac)
    T = Comrade.EHTLogClosureAmplitudeDatum{pol, eltype(lcamp.T), typeof(arrayconfig(dvis)[1])}
    lcamp_vals = Comrade.logclosure_amplitudes(measurement(dvis), cldmat)
    lcamp_sig = abs2.(noise(dvis) ./ measurement(dvis))
    lcamp_cov = cldmat * Diagonal(lcamp_sig) * transpose(cldmat)
    return Comrade.EHTObservationTable{T}(lcamp_vals, lcamp_cov, clac)
end


end # module
