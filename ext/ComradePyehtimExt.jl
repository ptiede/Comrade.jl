module ComradePyehtimExt

using Comrade
if isdefined(Base, :get_extension)
    using Pyehtim
    using StructArrays: StructVector, StructArray, append!!
    using LinearAlgebra
    using StaticArraysCore
else
    using ..Pyehtim
    using ..StructArrays: StructVector, StructArray, append!!
    using ..LinearAlgebra
    using ..StaticArraysCore
end


function getvisfield(obs)
    obsamps = obs.data
    u   = pyconvert(Vector,         obsamps["u"])
    v   = pyconvert(Vector,         obsamps["v"])
    err = pyconvert(Vector,         obsamps["sigma"])
    vis = pyconvert(Vector,         obsamps["vis"])
    t1  = pyconvert(Vector{Symbol}, obsamps["t1"])
    t2  = pyconvert(Vector{Symbol}, obsamps["t2"])
    time= pyconvert(Vector,         obsamps["time"])
    freq= fill(pyconvert(eltype(u), obs.rf), length(time))
    bw  = fill(pyconvert(eltype(u), obs.bw), length(time))
    baseline = tuple.(t1, t2)

    return  StructArray{Comrade.EHTVisibilityDatum{eltype(u)}}(
        measurement = vis,
        U = u,
        V = v,
        error = err,
        T = time,
        F = freq,
        bandwidth = bw,
        baseline = baseline
    )
end


function getampfield(obs)
    obsamps = obs.amp
    uamp   = pyconvert(Vector, obsamps["u"])
    vamp   = pyconvert(Vector, obsamps["v"])
    erramp = pyconvert(Vector, obsamps["sigma"])
    amps   = pyconvert(Vector, obsamps["amp"])
    t1     = pyconvert(Vector{Symbol}, obsamps["t1"])
    t2     = pyconvert(Vector{Symbol}, obsamps["t2"])
    baseline = tuple.(t1, t2)
    time = pyconvert(Vector, obsamps["time"])
    freq = fill(get_rf(obs), length(time))

    return  StructArray{Comrade.EHTVisibilityAmplitudeDatum{Float64}}(
        measurement = amps,
        U = uamp,
        V = vamp,
        error = erramp,
        T = time,
        F = freq,
        baseline = baseline
    )
end

function getcpfield(obs)
    obscp = obs.cphase
    u1 = pyconvert(Vector, obscp["u1"])
    v1 = pyconvert(Vector, obscp["v1"])
    u2 = pyconvert(Vector, obscp["u2"])
    v2 = pyconvert(Vector, obscp["v2"])
    u3 = pyconvert(Vector, obscp["u3"])
    v3 = pyconvert(Vector, obscp["v3"])
    cp    = deg2rad.(pyconvert(Vector, obscp["cphase"]))
    errcp = deg2rad.(pyconvert(Vector, obscp["sigmacp"]))

    t1 = pyconvert(Vector{Symbol}, obscp["t1"])
    t2 = pyconvert(Vector{Symbol}, obscp["t2"])
    t3 = pyconvert(Vector{Symbol}, obscp["t3"])
    baseline = tuple.(t1, t2, t3)
    time = pyconvert(Vector, obscp["time"])
    freq = fill(get_rf(obs), length(time))

    return StructArray{Comrade.EHTClosurePhaseDatum{Float64}}(
        measurement = cp,
        U1 = u1,
        V1 = v1,
        U2 = u2,
        V2 = v2,
        U3 = u3,
        V3 = v3,
        error = errcp,
        T = time,
        F = freq,
        triangle = baseline
    )

end

function getlcampfield(obs)
    obslcamp = obs.logcamp
    u1 = pyconvert(Vector, obslcamp["u1"])
    v1 = pyconvert(Vector, obslcamp["v1"])
    u2 = pyconvert(Vector, obslcamp["u2"])
    v2 = pyconvert(Vector, obslcamp["v2"])
    u3 = pyconvert(Vector, obslcamp["u3"])
    v3 = pyconvert(Vector, obslcamp["v3"])
    u4 = pyconvert(Vector, obslcamp["u4"])
    v4 = pyconvert(Vector, obslcamp["v4"])
    camp    = pyconvert(Vector, obslcamp["camp"])
    errcamp = pyconvert(Vector, obslcamp["sigmaca"])

    t1 = pyconvert(Vector{Symbol}, obslcamp["t1"])
    t2 = pyconvert(Vector{Symbol}, obslcamp["t2"])
    t3 = pyconvert(Vector{Symbol}, obslcamp["t3"])
    t4 = pyconvert(Vector{Symbol}, obslcamp["t4"])
    baseline = tuple.(t1, t2, t3, t4)
    time = pyconvert(Vector, obslcamp["time"])
    freq = fill(get_rf(obs), length(time))

    return StructArray{Comrade.EHTLogClosureAmplitudeDatum{Float64}}(
        measurement = camp,
        U1 = u1,
        V1 = v1,
        U2 = u2,
        V2 = v2,
        U3 = u3,
        V3 = v3,
        U4 = u4,
        V4 = v4,
        error = errcamp,
        T = time,
        F = freq,
        quadrangle = baseline
    )

end


function getcoherency(obs)

    # ensure that the visibilities are represented in a circular basis
    obs = obs.switch_polrep("circ")

    # get (u,v) coordinates
    u = pyconvert(Vector, obs.data["u"])
    v = pyconvert(Vector, obs.data["v"])

    # get visibilities
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


    # get timestamps and frequencies
    time = pyconvert(Vector, obs.data["time"])
    freq = fill(get_rf(obs), length(time))

    # get baseline info
    t1 = pyconvert(Vector{Symbol}, obs.data["t1"])
    t2 = pyconvert(Vector{Symbol}, obs.data["t2"])
    baseline = tuple.(t1, t2)

    # provide the polarization basis info
    single_polbasis = (CirBasis(), CirBasis())
    polbasis = fill(single_polbasis,length(u))

    # prepare output
    output = StructArray{Comrade.EHTCoherencyDatum{eltype(u),
                         typeof(single_polbasis[1]),
                         typeof(single_polbasis[2]),
                         eltype(cohmat),
                         eltype(errmat)}}(
        measurement = cohmat,
        error = errmat,
        U = u,
        V = v,
        T = time,
        F = freq,
        baseline = baseline,
        polbasis = polbasis
        )

    return output

end




"""
    extract_amp(obs)
Extracts the visibility amplitudes from an ehtim observation object.

Any valid keyword arguments to `add_amp` in ehtim can be passed through extract_amp.

Returns an EHTObservation with visibility amplitude data
"""
function Comrade.extract_amp(obsc; kwargs...)
    obs = obsc.copy()
    obs.add_scans()
    obs.reorder_tarr_snr()
    obs.add_amp(;kwargs...)
    data = getampfield(obs)
    ra, dec = get_radec(obs)
    mjd = get_mjd(obs)
    source = get_source(obs)
    bw = get_bw(obs)
    angles = get_fr_angles(obs)
    tarr = Pyehtim.get_arraytable(obsc)
    scans = get_scantable(obsc)
    ac = Comrade._arrayconfig(data, angles, tarr, scans, bw)
    return Comrade.EHTObservation(data = data, mjd = mjd,
                   ra = ra, dec= dec,
                   config = ac,
                   bandwidth=bw,
                   source = source,
    )
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
    extract_vis(obs; kwargs...)
Extracts the complex visibilities from an ehtim observation object

This grabs the raw `data` object from the obs object. Any keyword arguments are ignored.

Returns an EHTObservation with complex visibility data
"""
function Comrade.extract_vis(obsc; kwargs...)
    obs = obsc.copy()
    obs.add_scans()
    obs.reorder_tarr_snr()

    data = getvisfield(obs)
    ra, dec = get_radec(obs)
    mjd = get_mjd(obs)
    source = get_source(obs)
    bw = get_bw(obs)
    rf = get_rf(obs)
    angles = get_fr_angles(obs)
    tarr = Pyehtim.get_arraytable(obsc)
    sc = pyconvert(Matrix, obs.scans)
    scans = get_scantable(obs)
    ac = Comrade._arrayconfig(data, angles, tarr, scans, bw)
    return Comrade.EHTObservation(
                   data = data, mjd = mjd,
                   config=ac,
                   ra = ra, dec= dec,
                   bandwidth=bw,
                   source = source,
    )
end

"""
    extract_coherency(obs; kwargs...)
Extracts the coherency matrix from an ehtim observation object

This grabs the raw `data` object from the obs object. Any keyword arguments are ignored.

Returns an EHTObservation with coherency matrix
"""
function Comrade.extract_coherency(obsc; kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()
    obs.add_scans()
    data = getcoherency(obs)
    ra, dec = get_radec(obs)
    mjd = get_mjd(obs)
    source = get_source(obs)
    bw = get_bw(obs)
    rf = get_rf(obs)
    angles = get_fr_angles(obs)
    tarr = Pyehtim.get_arraytable(obs)
    scans = get_scantable(obs)
    ac = Comrade._arrayconfig(data, angles, tarr, scans, bw)
    return Comrade.EHTObservation(data = data, mjd = mjd,
                   ra = ra, dec= dec,
                   config = ac,
                   bandwidth=bw,
                   source = source,
    )
end



function closure_designmat(scancp::Comrade.Scan{A,B,C}, scanvis) where {A,B,C<:StructArray{<:Comrade.EHTClosurePhaseDatum}}
    ant1, ant2, ant3 = baselines(scancp)
    antvis1, antvis2 = baselines(scanvis)
    design_mat = zeros(length(scancp), length(scanvis))
    #throw("here")
    # fill in each row of the design matrix
    for i in axes(design_mat, 2), j in axes(design_mat,1)
        a1 = ant1[j]
        a2 = ant2[j]
        a3 = ant3[j]
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

function closure_designmat(scanlca::Comrade.Scan{A,B,C}, scanvis) where {A,B,C<:StructArray{<:Comrade.EHTLogClosureAmplitudeDatum}}
    ant1, ant2, ant3, ant4 = baselines(scanlca)
    antvis1, antvis2 = baselines(scanvis)
    design_mat = zeros(length(scanlca), length(scanvis))
    # fill in each row of the design matrix
    for i in axes(design_mat, 2), j in axes(design_mat,1)
        a1 = ant1[j]
        a2 = ant2[j]
        a3 = ant3[j]
        a4 = ant4[j]

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

    # now create EHTObservation
    ra, dec = get_radec(obs)
    mjd = get_mjd(obs)
    source = get_source(obs)
    bw = get_bw(obs)
    # ac = arrayconfig(dvis)
    clac = Comrade.ClosureConfig(dvis, dmat)
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
    data = getcpfield(obs)
    ra, dec = get_radec(obs)
    mjd = get_mjd(obs)
    source = get_source(obs)
    bw = get_bw(obs)


    cphase = Comrade.EHTObservation(data = data, mjd = mjd,
                   config=nothing,
                   ra = ra, dec= dec,
                   bandwidth=bw,
                   source = source,
    )

    stcp = scantable(cphase)

    #Now make the vis obs
    dvis = Comrade.extract_vis(obsc)
    st = scantable(dvis)
    S = eltype(dvis[:error])

    dmat = Matrix{S}[]
    for i in 1:length(st)
        scanvis = st[i]
        ind = findfirst(==(scanvis.time), stcp.times)
        if isnothing(ind)
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
        dmatscan = closure_designmat(stcp[ind], st[i])
        push!(dmat, dmatscan)
    end

    if iszero(dmat[1])
        dmat[2] = hcat(zeros(S, size(dmat[2],1), size(dmat[1],2)), dmat[2])
        dmat = dmat[2:end]
    end

    # ac = arrayconfig(dvis)
    clac = Comrade.ClosureConfig(dvis, dmat)
    return  Comrade.EHTObservation(data = data, mjd = mjd,
                           config=clac,
                           ra = ra, dec= dec,
                           bandwidth=bw,
                           source = source
                          )
end

function _make_lcamp(obsc, count="max"; kwargs...)
    obs = obsc.copy()
    obs.add_scans()
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

    lcamp = _make_lcamp(obsc; count, kwargs...)

    stlca = scantable(lcamp)

    #Now make the vis obs
    dvis = Comrade.extract_vis(obsc)
    st = scantable(dvis)
    S = eltype(dvis[:U])

    dmat = Matrix{S}[]
    for i in 1:length(st)
        scanvis = st[i]
        ind = findfirst(==(scanvis.time), stlca.times)
        if isnothing(ind)
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
        dmatscan = closure_designmat(stlca[ind], st[i])
        push!(dmat, dmatscan)
    end

    if iszero(dmat[1])
        dmat[2] = hcat(zeros(S, size(dmat[2],1), size(dmat[1],2)), dmat[2])
        dmat = dmat[2:end]
    end


    # ac = arrayconfig(dvis)
    clac = Comrade.ClosureConfig(dvis, dmat)
    return  Comrade.EHTObservation(data = lcamp.data, mjd = lcamp.mjd,
                           config=clac,
                           ra = lcamp.ra, dec= lcamp.dec,
                           bandwidth=lcamp.bandwidth,
                           source = lcamp.source
                          )
end


function minimal_cphase(obsc; kwargs...)
    # compute a maximum set of closure phases
    obs = obsc.copy()
    obs.add_scans()
    # reorder to maximize the snr
    obs.reorder_tarr_snr()

    cphase = _ehtim_cphase(obs; count="max", kwargs...)
    stcp = scantable(cphase)

    #Now make the vis obs
    dvis = Comrade.extract_vis(obsc)
    st = scantable(dvis)

    minset, dmat = _minimal_closure(stcp, st)
    # now create EHTObservation
    ra, dec = get_radec(obs)
    mjd = get_mjd(obs)
    source = get_source(obs)
    bw = get_bw(obs)
    # ac = arrayconfig(dvis)
    clac = Comrade.ClosureConfig(dvis, dmat)
    return Comrade.EHTObservation(data = minset, mjd = mjd,
                          config=clac,
                          ra = ra, dec= dec,
                          bandwidth=bw,
                          source = source,
                        )
end




function _minimal_closure(stcl, st)

    # Determine the number of timestamps containing a closure triangle
    T = typeof(stcl[1][1])
    S = typeof(stcl[1][1].U1)
    # loop over all timestamps
    minset = StructVector{T}(undef, 0)
    dmat = Matrix{S}[]
    for i in 1:length(st.times)

        scanvis = st[i]
        # find closure scan that matches time
        ind = findfirst(==(scanvis.time), stcl.times)
        if isnothing(ind)
            inow = length(dmat)
            # I want to construct block diagonal matrices for efficiency but we need
            # to be careful with scans that can't form closures. This hack moves these
            # scans into the preceding scan but fills it with zeros
            if i > 1
                dmat[inow] = hcat(dmat[inow], zeros(S, size(dmat[inow],1), length(scanvis.scan)))
            else
                push!(dmat, zeros(S, 1, length(scanvis.scan)))
            end
            continue
        end

        # Get our cphase scan
        scancl = stcl[ind]

        # sort by closure SNR so we form a nice minimal set
        snr = inv.(scancl.scan.error)
        ind_snr = sortperm(snr)
        scancl = scancl[ind_snr]
        snr = snr[ind_snr]

        # now find the visibility scan that matches the closure one based on timestamps


        # initialize the design matrix
        design_mat = closure_designmat(scancl, scanvis)

        # determine the expected size of the minimal set
        # this is needed to make sure we aren't killing too many triangles
        nmin = rank(design_mat)

        # print some info
        #println("For timestamp $(stcl.times[i]):")

        # get the current stations
        #println("Observing stations are $(stations(scancl))")

        #println("Size of maximal set of closure products = $(length(scancl))")
        #println("Size of minimal set of closure products = $(nmin)")
        #println("...")

        ##########################################################
        # start of loop to recover minimal set
        ##########################################################
        minset_scan, dmat_min = minimal_closure_scan(scancl, scanvis, nmin)
        minset = append!!(minset, minset_scan.scan)
        push!(dmat, dmat_min)
    end

    # Now if the first scan can't form a closure we will have an extra row of zeros
    # kill this
    if iszero(dmat[1])
        dmat[2] = hcat(zeros(S, size(dmat[2],1), size(dmat[1],2)), dmat[2])
        dmat = dmat[2:end]
    end
    return minset, dmat
end


function minimal_closure_scan(scancl0, scanvis, nmin::Int)
    # make a mask to keep track of which clhases will stick around
    scancl = deepcopy(scancl0)
    keep = fill(true, length(scancl))

    # remember the original minimal set size
    nmin0 = nmin

    # perform the loop
    count = 1
    keep = fill(true, length(scancl))
    good = true
    while good
        # recreate the mask each time
        scankeep = scancl[keep]

        design_mat = closure_designmat(scankeep, scanvis)

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
    scankeep = scancl[keep]
    dmat = closure_designmat(scankeep, scanvis)
    if length(scankeep) != nmin
        @error "minimal set not found $(length(scankeep)) $(nmin)"
        throw("No minimal set found at time $(scancl.time)")
    end
    return scankeep, dmat
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
function Comrade.extract_cphase(obs; count="min-correct", cut_trivial=false, uvmin=0.1e9,  kwargs...)
    if count == "min-correct"
        return minimal_cphase(obs; cut_trivial, uvmin, kwargs...)
    else
        return _ehtim_cphase(obs; count=count, cut_trivial=cut_trivial, uvmin=uvmin, kwargs...)
    end
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
function Comrade.extract_lcamp(obs; count="min-correct", kwargs...)
    if count == "min-correct"
        return minimal_lcamp(obs; kwargs...)
    else
        return _ehtim_lcamp(obs; count=count, kwargs...)
    end
end


end
