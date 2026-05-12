module ComradeVLBIFilesExt

using Comrade
using Accessors
using VLBIFiles
using VLBIFiles: VLBIData, VLBI
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
    VLBIData.AntennaMountType.AltAzimuth => (1.0, 0.0, 0.0),
    VLBIData.AntennaMountType.Equatorial => (0.0, 0.0, 0.0),
    VLBIData.AntennaMountType.Orbiting => (0.0, 0.0, 0.0),
    VLBIData.AntennaMountType.XY => (1.0, 0.0, 0.0),
    VLBIData.AntennaMountType.NaismithR => (1.0, 1.0, 0.0),
    VLBIData.AntennaMountType.NaismithL => (1.0, -1.0, 0.0),
    VLBIData.AntennaMountType.ApertureArray => (0.0, 0.0, 0.0),
    VLBIData.AntennaMountType.Unknown => (1.0, 0.0, 0.0),
)


function _mount_coeffs(ant)
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


function _mjd(date_obs)
    return Dates.value(Date(date_obs) - Date(1858, 11, 17))
end


function _build_tarr(
        ar::VLBIFiles.AntArray
    )
    # collect unique antennas across (typically just one) ant_array
    seen = Set{Symbol}()
    sites = Symbol[]
    X = Float64[]; Y = Float64[]; Z = Float64[]
    SEFD1 = Float64[]; SEFD2 = Float64[]
    fr_par = Float64[]; fr_el = Float64[]; fr_off = Float64[]
    for (_, ant) in pairs(ar.antennas)
        ant.name in seen && continue
        push!(seen, ant.name)
        push!(sites, ant.name)
        push!(X, ant.xyz[1])
        push!(Y, ant.xyz[2])
        push!(Z, ant.xyz[3])
        push!(SEFD1, 0.0)
        push!(SEFD2, 0.0)
        fp, fe, fo = _mount_coeffs(ant)
        push!(fr_par, fp); push!(fr_el, fe); push!(fr_off, fo)
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

function Comrade._preptable(obs::VLBIFiles.UVData, dataproduct::Comrade.VLBIDataProducts)

    time_average = get(dataproduct.keywords, :time_average, nothing)
    frequency_average = get(dataproduct.keywords, :frequency_average, true)

    uvtbl = VLBIData.uvtable(obs)
    antarray = only(obs.ant_arrays)
    ra, dec = _radec(obs)
    source = (; name = obs.header.object, ra = ra, dec = dec)
    puvtbl, scan_table = _prep_uvtable(
        uvtbl;
        time_average, frequency_average,
    )
    dataproduct2 = @set dataproduct.keywords = merge(NamedTuple(Comrade.keywords(dataproduct)), (; antarray, source, scan_table))
    return puvtbl, dataproduct2
end


function _build_scans(method, scan_table, uvtbl)
    if isnothing(scan_table)
        if hasproperty(uvtbl, :scan_id)
            sint = VLBIData.scan_intervals(uvtbl)
        else
            sint = VLBIData.scan_intervals(method, uvtbl)
        end
    else
        sint = scan_table
    end

    si = VLBIData.IntervalSets.endpoints.(sint)
    start = first.(si)
    stop = last.(si)
    day0 = Date(start[1])
    starthr = _hour_of_day.(start) .+ 24 .* Dates.value.(Date.(start) .- day0)
    stophr = _hour_of_day.(stop) .+ 24 .* Dates.value.(Date.(stop) .- day0)
    return StructArray((; start = starthr, stop = stophr))
end


function _arrayconfig(
        ptbl,
        uvtbl;
        antarray,
        source = (; name = "UNKNOWN", ra = 0.0, dec = 0.0),
        scan_table = nothing,
    )
    ra_deg, dec_deg = source.ra, source.dec
    # Pyehtim/ehtim stores RA in hours; match that convention so downstream
    # tooling that consumes either path agrees on units.
    ra_hours = ra_deg / 15
    mjd = _mjd(minimum(uvtbl.datetime))
    nm = source.name
    bw = _to_hz(uvtbl[1].freq_spec.width)
    tarr = _build_tarr(antarray)
    scans = _build_scans(VLBIData.GapBasedScans(), scan_table, uvtbl)

    n = length(ptbl)
    U = Vector{Float64}(undef, n)
    V = Vector{Float64}(undef, n)
    Ti = Vector{Float64}(undef, n)
    Fr = Vector{Float64}(undef, n)
    sites_v = Vector{Tuple{Symbol, Symbol}}(undef, n)
    elevation = StructVector{Tuple{Float64, Float64}}(undef, n)
    parallactic = StructVector{Tuple{Float64, Float64}}(undef, n)

    # caches built from antenna metadata
    site_xyz = Dict{Symbol, NTuple{3, Float64}}()
    poltype_map = Dict{Symbol, NTuple{2, Symbol}}()
    for (_, ant) in pairs(antarray.antennas)
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
    day0 = Date(minimum(ptbl.datetime))
    @inbounds for i in eachindex(ptbl)
        r = ptbl[i]
        s1, s2 = r.spec.bl.antennas
        U[i] = Float64(r.spec.uv.u)
        V[i] = Float64(r.spec.uv.v)
        Ti[i] = _hour_of_day(r.datetime) + 24 * Dates.value(Date(r.datetime) - day0)
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
        Symbol(nm), :UTC, data
    )
end


# polarization-row filtering helpers -------------------------------------

# Map Comrade pol Symbol -> the set of uvfits stokes labels and the linear
# combination weights to apply to (RR, LL, RL, LR) or (XX, YY, XY, YX) rows.
# Returns nothing if `pol` is one of the raw stokes labels (use directly).
function _pol_filter(pol::Symbol, available::Set{Symbol})
    if pol in available || pol == :all
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
        uvtbl;
        time_average = nothing, frequency_average = true
    )
    scan_table = VLBIData.scan_intervals(VLBI.GapBasedScans(), uvtbl)
    rows = uvtbl
    if frequency_average !== nothing && frequency_average !== false
        rows = VLBI.average_data(VLBI.ByFrequency(), rows)
    end
    if time_average !== nothing && time_average !== false
        rows = VLBI.average_data(time_average, rows)
    end

    return rows, scan_table
end

getpol(::Comrade.VLBIDataProducts) = PolarizedTypes.IPol
getpol(::Coherencies) = CoherencyMatrix


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
        uvtbl::AbstractArray{<:NamedTuple};
        antarray, source = (; name = "UNKNOWN", ra = 0.0, dec = 0.0),
        scan_table = nothing,
        pol = :I,
        kwargs...
    )


    val = VLBIData.uvtable_values_to(IPol, uvtbl)
    config = _arrayconfig(val, uvtbl; antarray, source, scan_table)
    vis, err = getvisfield(val)
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
        uvtbl::AbstractArray{<:NamedTuple};
        antarray, source = (; name = "UNKNOWN", ra = 0.0, dec = 0.0),
        scan_table = nothing,
        pol = :I,
        debias = false,
        kwargs...
    )

    val = VLBIData.uvtable_values_to(IPol, uvtbl)
    config = _arrayconfig(val, uvtbl; antarray, source, scan_table)
    vis, err = getvisfield(val)
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
        uvtbl::AbstractArray{<:NamedTuple};
        antarray, source = (; name = "UNKNOWN", ra = 0.0, dec = 0.0),
        scan_table = nothing,
        kwargs...
    )

    coherency = VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl)
    config = _arrayconfig(coherency, uvtbl; antarray, source, scan_table)

    cohmat = StructArray{SMatrix{2, 2, ComplexF64, 4}}(undef, length(coherency))
    errmat = StructArray{SMatrix{2, 2, Float64, 4}}(undef, length(coherency))

    for i in eachindex(coherency)
        c = coherency[i].value
        v = getproperty.(c, :v)
        e = getproperty.(c, :u)
        cohmat[i] = v
        errmat[i] = e
    end

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
        uvtbl::AbstractArray{<:NamedTuple};
        antarray, source = (; name = "UNKNOWN", ra = 0.0, dec = 0.0),
        pol = :I, count = "min",
        kwargs...
    )
    dvis = Comrade.extract_vis(
        uvtbl; antarray, source, pol,
        kwargs...
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
        uvtbl::AbstractArray{<:NamedTuple};
        antarray, source = (; name = "UNKNOWN", ra = 0.0, dec = 0.0),
        pol = :I, count = "min",
        kwargs...
    )
    dvis = Comrade.extract_vis(
        uvtbl; antarray, source, pol,
        kwargs...
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
