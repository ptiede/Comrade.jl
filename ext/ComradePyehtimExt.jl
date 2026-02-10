module ComradePyehtimExt

using Comrade
using Pyehtim
using StructArrays: StructVector, StructArray, append!!
using LinearAlgebra
using StaticArraysCore


removeesc(str) = Symbol(replace(str, r"[^a-zA-Z0-9\s]" => ""))

function build_arrayconfig(obsc)
    obsd = obsc.data
    ra, dec = get_radec(obsc)
    mjd = get_mjd(obsc)
    source = get_source(obsc)
    bw = get_bw(obsc)
    angles = get_fr_angles(obsc)
    tarr = StructArray(Pyehtim.get_arraytable(obsc))

    # This is because sometimes eht-imaging sets the scans to nothing
    # and sometimes it fills it with junk
    PythonCall.Core.pyisnone(obsc.scans) && obsc.add_scans()
    if length(obsc.scans) <= 1
        obsc.add_scans()
    end


    scans = StructArray(get_scantable(obsc))
    bw = get_bw(obsc)
    elevation = StructArray(angles[1])
    parallactic = StructArray(angles[2])

    U = pyconvert(Vector, obsd["u"])
    V = pyconvert(Vector, obsd["v"])
    t1 = removeesc.(pyconvert(Vector{String}, obsd["t1"]))
    t2 = removeesc.(pyconvert(Vector{String}, obsd["t2"]))
    Ti = pyconvert(Vector, obsd["time"])
    Fr = fill(pyconvert(eltype(U), obsc.rf), length(U))
    sites = tuple.(t1, t2)
    single_polbasis = (CirBasis(), CirBasis())
    polbasis = fill(single_polbasis, length(U))
    data = StructArray{Comrade.EHTArrayBaselineDatum{eltype(U), eltype(polbasis), eltype(elevation[1][1])}}(
        (; U, V, Ti, Fr, sites, polbasis, elevation, parallactic)
    )
    return Comrade.EHTArrayConfiguration(bw, tarr, scans, mjd, ra, dec, source, :UTC, data)
end


function getvisfield(obs)
    obsd = obs.data
    vis = pyconvert(Vector{ComplexF64}, obsd["vis"])
    err = pyconvert(Vector{Float64}, obsd["sigma"])
    return vis, err
end


function getampfield(obs)
    obsamps = obs.amp
    erramp = pyconvert(Vector, obsamps["sigma"])
    amps = pyconvert(Vector, obsamps["amp"])
    return amps, erramp
end

function getcoherency(obs)

    # check if the obs is in circular basis otherwise noise out
    @assert(
        (pyconvert(String, obs.polrep) == "circ"),
        "obs is not in circular polarization.\n" *
            "To fix read in the data using\n" *
            "  Pyehtim.load_uvfits_and_array(obsname, arrayname, polrep=\"circ\")\n" *
            "Do not use\n  obs.switch_polrep(\"circ\")\nsince missing hands will not be handled correctly."
    )


    c11 = pyconvert(Vector, obs.data["rrvis"])
    c12 = pyconvert(Vector, obs.data["rlvis"])
    c21 = pyconvert(Vector, obs.data["lrvis"])
    c22 = pyconvert(Vector, obs.data["llvis"])

    cohmat = StructArray{SMatrix{2, 2, eltype(c11), 4}}((c11, c21, c12, c22))

    # get uncertainties
    e11 = copy(pyconvert(Vector, obs.data["rrsigma"]))
    e12 = copy(pyconvert(Vector, obs.data["rlsigma"]))
    e21 = copy(pyconvert(Vector, obs.data["lrsigma"]))
    e22 = copy(pyconvert(Vector, obs.data["llsigma"]))

    errmat = StructArray{SMatrix{2, 2, eltype(e11), 4}}((e11, e21, e12, e22))

    return cohmat, errmat

end


function getcpfield(obs)
    # Here we just return the information needed to form
    # the closure configuration
    obscp = obs.cphase
    time = pyconvert(Vector, obscp["time"])
    freq = fill(get_rf(obs), length(time))
    t1 = removeesc(pyconvert(Vector{String}, obscp["t1"]))
    t2 = removeesc(pyconvert(Vector{String}, obscp["t2"]))
    t3 = removeesc(pyconvert(Vector{String}, obscp["t3"]))
    noise = pyconvert(Vector, obscp["sigmacp"])
    baseline = tuple.(t1, t2, t3)
    return StructArray((; T = time, F = freq, noise, baseline))
end

function getlcampfield(obs)
    # Here we just return the information needed to form
    # the closure configuration
    obslcamp = obs.logcamp
    t1 = removeesc(pyconvert(Vector{String}, obslcamp["t1"]))
    t2 = removeesc(pyconvert(Vector{String}, obslcamp["t2"]))
    t3 = removeesc(pyconvert(Vector{String}, obslcamp["t3"]))
    time = pyconvert(Vector, obslcamp["time"])
    t4 = removeesc(pyconvert(Vector{String}, obslcamp["t4"]))
    baseline = tuple.(t1, t2, t3, t4)
    noise = pyconvert(Vector, obslcamp["sigmaca"])
    freq = fill(get_rf(obs), length(time))
    return StructArray((; T = time, F = freq, baseline, noise))
end


function get_arraytable(obs)
    return StructArray(
        sites = removeesc(pyconvert(Vector{String}, obs.tarr["site"])),
        X = pyconvert(Vector, obs.tarr["x"]),
        Y = pyconvert(Vector, obs.tarr["y"]),
        Z = pyconvert(Vector, obs.tarr["z"]),
        SEFD1 = pyconvert(Vector, obs.tarr["sefdr"]),
        SEFD2 = pyconvert(Vector, obs.tarr["sefdl"]),
        fr_parallactic = pyconvert(Vector, obs.tarr["fr_par"]),
        fr_elevation = pyconvert(Vector, obs.tarr["fr_elev"]),
        fr_offset = deg2rad.(pyconvert(Vector, obs.tarr["fr_off"])),
    )
end


"""
    extract_amp(obs)
Extracts the visibility amplitudes from an ehtim observation object.

Any valid keyword arguments to `add_amp` in ehtim can be passed through extract_amp.

Returns an EHTObservationTable with visibility amplitude data
"""
function Comrade.extract_amp(obsc; pol = :I, debias = false, kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()
    obs.add_amp(; debias, kwargs...)
    config = build_arrayconfig(obs)
    amp, amperr = getampfield(obs)
    T = Comrade.EHTVisibilityAmplitudeDatum{pol, eltype(amp), typeof(config[1])}
    return Comrade.EHTObservationTable{T}(amp, amperr, config)
end


"""
    extract_vis(obs; kwargs...)
Extracts the complex visibilities from an ehtim observation object

This grabs the raw `data` object from the obs object. Any keyword arguments are ignored.

Returns an EHTObservationTable with complex visibility data
"""
function Comrade.extract_vis(obsc; pol = :I, kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()
    config = build_arrayconfig(obs)
    vis, viserr = getvisfield(obs)
    T = Comrade.EHTVisibilityDatum{pol, eltype(viserr), typeof(config[1])}
    return Comrade.EHTObservationTable{T}(vis, viserr, config)
end

"""
    extract_coherency(obs; kwargs...)
Extracts the coherency matrix from an ehtim observation object

This grabs the raw `data` object from the obs object. Any keyword arguments are ignored.

Returns an EHTObservationTable with coherency matrices
"""
function Comrade.extract_coherency(obsc; kwargs...)
    obs = obsc.copy()
    config = build_arrayconfig(obs)
    vis, viserr = getcoherency(obs)
    T = Comrade.EHTCoherencyDatum{eltype(real(vis[1])), typeof(config[1]), eltype(vis), eltype(viserr)}
    return Comrade.EHTObservationTable{T}(vis, viserr, config)
end

function _ehtim_cphase(obsc; count = "max", cut_trivial = false, uvmin = 0.1e9, kwargs...)
    obs = obsc.copy()

    # cut 0 baselines since these are trivial triangles
    if cut_trivial
        obs = obs.flag_uvdist(uv_min = uvmin)
    end

    obs.reorder_tarr_snr()

    obs.add_cphase(; count = count, kwargs...)
    cphase = getcpfield(obs)
    return cphase, Comrade.extract_vis(obs)
end

function _make_lcamp(obsc, count = "max"; kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()

    obs.add_logcamp(; count = count, kwargs...)
    data = getlcampfield(obs)
    ra, dec = get_radec(obs)
    mjd = get_mjd(obs)
    source = get_source(obs)
    bw = get_bw(obs)

    return Comrade.EHTObservation(
        data = data, mjd = mjd,
        config = nothing,
        ra = ra, dec = dec,
        bandwidth = bw,
        source = source,
    )

end


function _ehtim_lcamp(obsc; count = "max", kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()

    obs.add_logcamp(; count = count, kwargs...)
    lcamp = getlcampfield(obs)
    return lcamp, Comrade.extract_vis(obs)
end


"""
    extract_cphase(obs)
Extracts the closure phases from an ehtim observation object

Any valid keyword arguments to `add_cphase` in ehtim can be passed through extract_cphase.

Returns an EHTObservation with closure phases datums

# Special Keyword arguments:
 - count: How the closures are formed, the available options are "min-correct", "min", "max"
 - cut_trivial: Cut the trivial triangles from the closures
 - uvmin: The flag to decide what are trivial triangles. Any baseline with ||(u,v)|| < uvmin
          are removed.
 - kwargs...: Other arguments are forwarded to eht-imaging.

# Warning

The `count` keyword argument is treated specially in `Comrade`. The default option
is "min-correct" and should almost always be used.
This option construct a minimal set of closure phases that is valid even when
the array isn't fully connected. For testing and legacy reasons we `ehtim` other count
options are also included. However, the current `ehtim` count="min" option is broken
and does construct proper minimal sets of closure quantities if the array isn't fully connected.

"""
function Comrade.extract_cphase(obs; pol = :I, count = "min", kwargs...)
    obsc = obs.copy()
    obsc.reorder_tarr_snr()

    cphase, dvis = _ehtim_cphase(obsc; count = "max", kwargs...)
    clac = Comrade.build_closure_config(dvis, cphase; type = :cphase, count)
    T = Comrade.EHTClosurePhaseDatum{pol, eltype(cphase.T), typeof(arrayconfig(dvis)[1])}
    cp = Comrade.closure_phases(measurement(dvis), Comrade.designmat(clac))
    cp_sig = abs2.(Comrade.noise(dvis) ./ Comrade.measurement(dvis))
    cp_cov = Comrade.designmat(clac) * Diagonal(cp_sig) * transpose(Comrade.designmat(clac))
    return Comrade.EHTObservationTable{T}(cp, cp_cov, clac)
end


"""
    extract_lcamp(obs; kwargs...)
Extracts the log-closure amp. from an ehtim observation object

Any valid keyword arguments to `add_logcamp` in ehtim can be passed through extract_lcamp.

# Special Keyword arguments:
 - count: How the closures are formed, the available options are "min-correct", "min", "max"
 - kwargs...: Other arguments are forwarded to eht-imaging.

Returns an EHTObservation with log-closure amp. datums

# Warning
The `count` keyword argument is treated specially in `Comrade`. The default option
is "min-correct" and should almost always be used.
This option construct a minimal set of closure phases that is valid even when
the array isn't fully connected. For testing and legacy reasons we `ehtim` other count
options are also included. However, the current `ehtim` count="min" option is broken
and does construct proper minimal sets of closure quantities if the array isn't fully connected.

"""
function Comrade.extract_lcamp(obs; pol = :I, count = "min", kwargs...)
    obsc = obs.copy()
    obsc.reorder_tarr_snr()

    lcamp, dvis = _ehtim_lcamp(obsc; count = "max", kwargs...)
    clac = Comrade.build_closure_config(dvis, lcamp; type = :lcamp, count)
    cldmat = Comrade.designmat(clac)
    T = Comrade.EHTLogClosureAmplitudeDatum{pol, eltype(lcamp.T), typeof(arrayconfig(dvis)[1])}
    lcamp_vals = Comrade.logclosure_amplitudes(measurement(dvis), cldmat)
    lcamp_sig = abs2.(Comrade.noise(dvis) ./ Comrade.measurement(dvis))
    lcamp_cov = cldmat * Diagonal(lcamp_sig) * transpose(cldmat)
    return Comrade.EHTObservationTable{T}(lcamp_vals, lcamp_cov, clac)
end


end
