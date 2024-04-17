module ComradePyehtimExt

using Comrade
if isdefined(Base, :get_extension)
    using Pyehtim
    using StructArrays: StructVector, StructArray, append!!
    using LinearAlgebra
    using StaticArraysCore
    using TypedTables
else
    using ..Pyehtim
    using ..StructArrays: StructVector, StructArray, append!!
    using ..LinearAlgebra
    using ..StaticArraysCore
    using ..TypedTables
end

function build_arrayconfig(obs)
    obsd = obs.data
    obsc = obs.copy()
    ra, dec = get_radec(obsc)
    mjd = get_mjd(obsc)
    source = get_source(obsc)
    bw = get_bw(obsc)
    angles = get_fr_angles(obsc)
    tarr = Pyehtim.get_arraytable(obsc)
    scans = get_scantable(obsc)
    bw  = get_bw(obsc)
    elevation = StructArray(angles[1])
    parallactic  = StructArray(angles[2])

    U   = pyconvert(Vector,         obsd["u"])
    V   = pyconvert(Vector,         obsd["v"])
    t1  = pyconvert(Vector{Symbol}, obsd["t1"])
    t2  = pyconvert(Vector{Symbol}, obsd["t2"])
    T= pyconvert(Vector,         obsd["time"])
    F= fill(pyconvert(eltype(U), obsc.rf), length(U))
    sites = tuple.(t1, t2)
    single_polbasis = (CirBasis(), CirBasis())
    polbasis = fill(single_polbasis,length(U))
    data = StructArray{Comrade.EHTArrayBaselineDatum{eltype(U), eltype(polbasis), eltype(elevation[1][1])}}(
                (;U, V, T, F, sites, polbasis, elevation, parallactic)
    )
    return Comrade.EHTArrayConfiguration(bw, tarr, scans, mjd, ra, dec, source, :UTC, data)
end


function getvisfield(obs)
    obsd = obs.data
    vis  = pyconvert(Vector{ComplexF64}, obsd["vis"])
    err  = pyconvert(Vector{Float64}, obsd["sigma"])
    return vis, err
end


function getampfield(obs)
    obsamps = obs.amp
    erramp = pyconvert(Vector, obsamps["sigma"])
    amps   = pyconvert(Vector, obsamps["amp"])
    return amps, erramp
end

function getcoherency(obs)

    # check if the obs is in circular basis otherwise noise out
    @assert((pyconvert(String, obs.polrep) == "circ"),
            "obs is not in circular polarization.\n"*
            "To fix read in the data using\n"*
            "  Pyehtim.load_uvfits_and_array(obsname, arrayname, polrep=\"circ\")\n"*
            "Do not use\n  obs.switch_polrep(\"circ\")\nsince missing hands will not be handled correctly."
        )


    c11 = pyconvert(Vector, obs.data["rrvis"])
    c12 = pyconvert(Vector, obs.data["rlvis"])
    c21 = pyconvert(Vector, obs.data["lrvis"])
    c22 = pyconvert(Vector, obs.data["llvis"])

    cohmat = StructArray{SMatrix{2,2,eltype(c11), 4}}((c11, c21, c12, c22))

    # get uncertainties
    e11 = copy(pyconvert(Vector, obs.data["rrsigma"]))
    e12 = copy(pyconvert(Vector, obs.data["rlsigma"]))
    e21 = copy(pyconvert(Vector, obs.data["lrsigma"]))
    e22 = copy(pyconvert(Vector, obs.data["llsigma"]))

    errmat = StructArray{SMatrix{2,2,eltype(e11), 4}}((e11, e21, e12, e22))

    return cohmat, errmat

end


function getcpfield(obs)
    # Here we just return the information needed to form
    # the closure configuration
    obscp = obs.cphase
    time = pyconvert(Vector, obscp["time"])
    freq = fill(get_rf(obs), length(time))
    t1 = pyconvert(Vector{Symbol}, obscp["t1"])
    t2 = pyconvert(Vector{Symbol}, obscp["t2"])
    t3 = pyconvert(Vector{Symbol}, obscp["t3"])
    noise = pyconvert(Vector, obscp["sigmacp"])
    baseline = tuple.(t1, t2, t3)
    return Table((;T=time, F=freq, noise, baseline))
end

function getlcampfield(obs)
    # Here we just return the information needed to form
    # the closure configuration
    obslcamp = obs.logcamp
    t1 = pyconvert(Vector{Symbol}, obslcamp["t1"])
    t2 = pyconvert(Vector{Symbol}, obslcamp["t2"])
    t3 = pyconvert(Vector{Symbol}, obslcamp["t3"])
    time = pyconvert(Vector, obslcamp["time"])
    t4 = pyconvert(Vector{Symbol}, obslcamp["t4"])
    baseline = tuple.(t1, t2, t3, t4)
    noise = pyconvert(Vector, obslcamp["sigmaca"])
    freq = fill(get_rf(obs), length(time))
    return Table((;T=time, F=freq, baseline, noise))
end


function get_arraytable(obs)
    return Table(
        sites = pyconvert(Vector{Symbol}, obs.tarr["site"]),
        X     = pyconvert(Vector, obs.tarr["x"]),
        Y     = pyconvert(Vector, obs.tarr["y"]),
        Z     = pyconvert(Vector, obs.tarr["z"]),
        SEFD1 = pyconvert(Vector, obs.tarr["sefdr"]),
        SEFD2 = pyconvert(Vector, obs.tarr["sefdl"]),
        fr_parallactic = pyconvert(Vector, obs.tarr["fr_par"]),
        fr_elevation   = pyconvert(Vector, obs.tarr["fr_elev"]),
        fr_offset      = deg2rad.(pyconvert(Vector, obs.tarr["fr_off"])),
    )
end





"""
    extract_amp(obs)
Extracts the visibility amplitudes from an ehtim observation object.

Any valid keyword arguments to `add_amp` in ehtim can be passed through extract_amp.

Returns an EHTObservationTable with visibility amplitude data
"""
function Comrade.extract_amp(obsc; kwargs...)
    obs = obsc.copy()
    obs.add_scans()
    obs.reorder_tarr_snr()
    obs.add_amp(;kwargs...)
    config = build_arrayconfig(obs)
    amp, amperr = getampfield(obs)
    T = Comrade.EHTVisibilityAmplitudeDatum{eltype(amp), typeof(config[1])}
    return Comrade.EHTObservationTable{T}(amp, amperr, config)
end





"""
    extract_vis(obs; kwargs...)
Extracts the complex visibilities from an ehtim observation object

This grabs the raw `data` object from the obs object. Any keyword arguments are ignored.

Returns an EHTObservationTable with complex visibility data
"""
function Comrade.extract_vis(obsc; kwargs...)
    obs = obsc.copy()
    obs.add_scans()
    obs.reorder_tarr_snr()
    config = build_arrayconfig(obs)
    vis, viserr = getvisfield(obs)
    T = Comrade.EHTComplexVisibilityDatum{eltype(viserr), typeof(config[1])}
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
    obs.add_scans()
    obs.reorder_tarr_snr()
    config = build_arrayconfig(obs)
    vis, viserr = getcoherency(obs)
    T = Comrade.EHTCoherencyDatum{eltype(real(vis[1])), typeof(config[1]), eltype(vis), eltype(viserr)}
    return Comrade.EHTObservationTable{T}(vis, viserr, config)
end

function closure_designmat(type, closures, scanvis)
    if type == :cphase
        design_mat = closurephase_designmat(closures, scanvis)
    elseif type == :lcamp
        design_mat = closureamp_designmat(closures, scanvis)
    else
        throw(ArgumentError("Not a valid type of closure"))
    end
end


function closurephase_designmat(cphase, scanvis)
    antvis1, antvis2 = baselines(scanvis)
    design_mat = zeros(length(cphase), length(scanvis))
    #throw("here")
    # fill in each row of the design matrix
    for i in axes(design_mat, 2), j in axes(design_mat,1)
        a1, a2, a3 = cphase[j].baseline
        # check leg 1
        ((antvis1[i] == a1) & (antvis2[i] == a2)) && (design_mat[j,i] = 1.0)
        ((antvis1[i] == a2) & (antvis2[i] == a1)) && (design_mat[j,i] = -1.0)

        #check leg 2
        ((antvis1[i] == a2) & (antvis2[i] == a3)) && (design_mat[j,i] = 1.0)
        ((antvis1[i] == a3) & (antvis2[i] == a2)) && (design_mat[j,i] = -1.0)

        #check leg 3
        ((antvis1[i] == a3) & (antvis2[i] == a1)) && (design_mat[j,i] = 1.0)
        ((antvis1[i] == a1) & (antvis2[i] == a3)) && (design_mat[j,i] = -1.0)
    end
    return design_mat
end

function closureamp_designmat(lcamp, scanvis)
    antvis1, antvis2 = baselines(scanvis)
    design_mat = zeros(length(lcamp), length(scanvis))
    # fill in each row of the design matrix
    for i in axes(design_mat, 2), j in axes(design_mat,1)
        a1, a2, a3, a4 = lcamp[j].baseline

        av1 = antvis1[i]
        av2 = antvis2[i]
        # check leg 1
        (((av1 == a1)&(av2 == a2)) || (av1 == a2)&(av2 == a1))&& (design_mat[j,i] = 1.0)

        #check leg 2
        (((av1 == a3)&(av2 == a4)) || (av1 == a4)&(av2 == a3))&& (design_mat[j,i] = 1.0)

        #check leg 3
        (((av1 == a1)&(av2 == a4)) || (av1 == a4)&(av2 == a1))&& (design_mat[j,i] = -1.0)

        #check leg 4
        (((av1 == a2)&(av2 == a3)) || (av1 == a3)&(av2 == a2))&& (design_mat[j,i] = -1.0)

    end
    return design_mat
end

function build_dmats(type::Symbol, closure, st)
    S = eltype(closure.time)
    dmat = Matrix{S}[]
    for i in 1:length(st)
        scanvis = st[i]
        inds = findall(==(scanvis.time), closure.times)
        if isnothing(inds)
            inow = length(dmat)
            # I want to construct block diagonal matrices for efficienty but we need
            # to be careful with scan that can't form closures. This hack moves these
            # scans into the preceding scan but fills it with zeros
            if i > 1
                dmat[inow] = hcat(dmat[inow], zeros(S, size(dmat[inow],1), length(scanvis.scan)))
            else
                push!(dmat, zeros(S, 1, length(scanvis.scan)))
            end
            continue
        end
        scancl = closures[inds]
        if type == :cphase
            dmatscan = closurephase_designmat(scancl, scanvis)
        elseif type == :lcamp
            dmatscan = closureamp_designmat(scancl, scanvis)
        else
            throw(ArgumentError("Not a valid type of closure"))
        end
        push!(dmat, dmatscan)
    end

    if iszero(dmat[1])
        dmat[2] = hcat(zeros(S, size(dmat[2],1), size(dmat[1],2)), dmat[2])
        dmat = dmat[2:end]
    end
    return dmat
end


function minimal_lcamp(obsc; kwargs...)

    obs = obsc.copy()
    obs.add_scans()
    # reorder to maximize the snr
    obs.reorder_tarr_snr()

    lcamp = _ehtim_lcamp(obs; count="max", kwargs...)
    stlca = scantable(lcamp)

    #Now make the vis obs
    dvis = Comrade.extract_vis(obsc)
    st = scantable(dvis)

    minset, dmat = _minimal_closure(stlca, st)

    config = build_arrayconfig(obs)
    clac = Comrade.ClosureConfig(config, dmat)
    return Comrade.EHTObservation(data = minset, mjd = mjd,
                          config=clac,
                          ra = ra, dec= dec,
                          bandwidth=bw,
                          source = source,
                        )

end


function _ehtim_cphase(obsc; count="max", cut_trivial=false, uvmin=0.1e9, kwargs...)
    obs = obsc.copy()

    # cut 0 baselines since these are trivial triangles
    if cut_trivial
        obs = obs.flag_uvdist(uv_min=uvmin)
    end

    obs.reorder_tarr_snr()

    obs.add_cphase(;count=count, kwargs...)
    cphase = getcpfield(obs)
    return cphase, Comrade.extract_vis(obs)
end

function _make_lcamp(obsc, count="max"; kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()

    obs.add_logcamp(;count=count, kwargs...)
    data = getlcampfield(obs)
    ra, dec = get_radec(obs)
    mjd = get_mjd(obs)
    source = get_source(obs)
    bw = get_bw(obs)

    return Comrade.EHTObservation(data = data, mjd = mjd,
                   config=nothing,
                   ra = ra, dec= dec,
                   bandwidth=bw,
                   source = source,
    )

end


function _ehtim_lcamp(obsc; count="max", kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()

    obs.add_logcamp(;count=count, kwargs...)
    lcamp = getlcampfield(obs)
    return lcamp, Comrade.extract_vis(obs)
end



function _minimal_closure(type, closures, st)

    # Determine the number of timestamps containing a closure triangle
    S = eltype(closures.T)
    # loop over all timestamps
    dmat = Matrix{S}[]
    for i in 1:length(st)
        scanvis = st[i]
        inds = findall(==(scanvis.time), closures.T)
        if length(inds) == 0
            inow = length(dmat)
            # I want to construct block diagonal matrices for efficienty but we need
            # to be careful with scan that can't form closures. This hack moves these
            # scans into the preceding scan but fills it with zeros
            if i > 1
                dmat[inow] = hcat(dmat[inow], zeros(S, size(dmat[inow],1), length(scanvis.scan)))
            else
                push!(dmat, zeros(S, 1, length(scanvis.scan)))
            end
            continue
        end

        # Get our cphase scan
        scancl = closures[inds]
        # @info scancl

        # sort by closure noise so we form a nice minimal set
        snr = inv.(scancl.noise)
        ind_snr = sortperm(snr)
        scancl = scancl[ind_snr]

        # initialize the design matrix
        design_mat = closure_designmat(type, scancl, scanvis)
        # determine the expected size of the minimal set
        # this is needed to make sure we aren't killing too many triangles
        nmin = rank(design_mat)

        # print some info
        # println("For timestamp $(scanvis.time):")

        # get the current sites
        # println("Observing sites are $(scancl.baseline)")
        # println("scanvis sites ", sites(scanvis))

        # println("Size of maximal set of closure products = $(length(scancl))")
        # println("Size of minimal set of closure products = $(nmin)")
        # println("...")

        ##########################################################
        # start of loop to recover minimal set
        ##########################################################
        dmat_min = minimal_closure_scan(type, scancl, scanvis, nmin)
        push!(dmat, dmat_min)
    end

    # Now if the first scan can't form a closure we will have an extra row of zeros
    # kill this
    if iszero(dmat[1])
        dmat[2] = hcat(zeros(S, size(dmat[2],1), size(dmat[1],2)), dmat[2])
        dmat = dmat[2:end]
    end
    return dmat
end


function minimal_closure_scan(type, closures, scanvis, nmin::Int)
    # make a mask to keep track of which clhases will stick around
    scancl = deepcopy(closures)
    keep = fill(true, length(closures))

    # remember the original minimal set size
    nmin0 = nmin

    # perform the loop
    count = 1
    keep = fill(true, length(closures))
    good = true
    while good
        # recreate the mask each time
        closurekeep = closures[keep]

        design_mat = closure_designmat(type, closurekeep, scanvis)

        # determine the size of the minimal set
        nmin = rank(design_mat)
        #println(nmin0, " ", nmin, " ", size(design_mat, 1))

        # Now decide whether to continue pruning
        if (sum(keep) == nmin0) && (nmin == nmin0)
            good = false
        else
            if nmin == nmin0
                keep[count] = false
            else
                keep[count-1] = true
                count -= 1
            end
            #(nmin == nmin0) && (keep[count] = false)
        end
        count += 1
        count+1 > length(scancl) && break
    end

    # print out the size of the recovered set for double-checking
    closurekeep = closures[keep]
    dmat = closure_designmat(type, closurekeep, scanvis)
    if length(closurekeep) != nmin
        @error "minimal set not found $(length(closurekeep)) $(nmin)"
        throw("No minimal set found at time $(scanvis.time)")
    end
    return dmat
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
function Comrade.extract_cphase(obs; count="min", kwargs...)
    # compute a maximum set of closure phases
    obsc = obs.copy()
    # reorder to maximize the snr
    obsc.reorder_tarr_snr()

    cphase, dvis = _ehtim_cphase(obsc; count="max", kwargs...)

    #Now make the vis obs
    st = scantable(dvis)

    if count == "min"
        dmat = _minimal_closure(:cphase, cphase, st)
    elseif count == "max"
        dmat = build_dmats(:cphase, cphase, st)
    else
        throw(ArgumentError("$(count) is not valid use 'min' or 'max'"))
    end
    clac = Comrade.ClosureConfig(arrayconfig(dvis), dmat, measurement(dvis), noise(dvis))
    T = Comrade.EHTClosurePhaseDatum{eltype(cphase.T), typeof(arrayconfig(dvis)[1])}
    cp = Comrade.closure_phases(measurement(dvis), clac)
    cp_sig = abs2.(Comrade.noise(dvis)./Comrade.measurement(dvis))
    cp_cov = Comrade.designmat(clac)*Diagonal(cp_sig)*transpose(Comrade.designmat(clac))
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
function Comrade.extract_lcamp(obs; count="min", kwargs...)
    # compute a maximum set of closure phases
    obsc = obs.copy()
    # reorder to maximize the snr
    obsc.reorder_tarr_snr()

    lcamp, dvis = _ehtim_lcamp(obsc; count="max", kwargs...)

    #Now make the vis obs
    st = scantable(dvis)

    if count == "min"
        dmat = _minimal_closure(:lcamp, lcamp, st)
    elseif count == "max"
        dmat = build_dmats(:lcamp, lcamp, st)
    else
        throw(ArgumentError("$(count) is not valid use 'min' or 'max'"))
    end
    clac = Comrade.ClosureConfig(arrayconfig(dvis), dmat, measurement(dvis), noise(dvis))
    cldmat = Comrade.designmat(clac)
    T = Comrade.EHTLogClosureAmplitudeDatum{eltype(lcamp.T), typeof(arrayconfig(dvis)[1])}
    lcamp = Comrade.logclosure_amplitudes(measurement(dvis), clac)
    lcamp_sig = abs2.(Comrade.noise(dvis)./Comrade.measurement(dvis))
    lcamp_cov = cldmat*Diagonal(lcamp_sig)*transpose(cldmat)
    return Comrade.EHTObservationTable{T}(lcamp, lcamp_cov, clac)
end


end
