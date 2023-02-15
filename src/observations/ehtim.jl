export extract_amp, extract_vis, extract_lcamp, extract_cphase,
       extract_coherency,
       load_ehtim_uvfits, scan_average
using PyCall: set!

"""
    load_ehtim_uvfits(uvfile, arrayfile=nothing; kwargs...)

Load a uvfits file with eht-imaging and returns a eht-imaging `Obsdata`
object. You can optionally pass an array file as well that will load
additional information such at the telescopes field rotation information
with the arrayfile. This is expected to be an eht-imaging produced array
or antenna file.
"""
function load_ehtim_uvfits(uvfile, arrayfile=nothing; kwargs...)
    ehtim == PyNULL() && load_ehtim()
    obs = ehtim.obsdata.load_uvfits(uvfile; kwargs...)
    if arrayfile !== nothing
        tarr = ehtim.io.load.load_array_txt(arrayfile).tarr
        obs.tarr = tarr
    end
    return obs
end

function getvisfield(obs)
    obsamps = obs.data::PyObject
    u = deepcopy((get(obsamps, Vector{Float64}, "u")))
    v = deepcopy((get(obsamps, Vector{Float64}, "v")))
    err = deepcopy((get(obsamps, Vector{Float64}, "sigma")))
    vis = deepcopy((get(obsamps, Vector{Complex{Float64}}, "vis")))
    t1 = Symbol.(deepcopy((get(obsamps, Vector{String}, "t1"))))
    t2 = Symbol.(deepcopy((get(obsamps, Vector{String}, "t2"))))
    baseline = tuple.(t1, t2)
    time = deepcopy((get(obsamps, Vector{Float64}, "time")))
    freq = fill(obs.rf, length(time))
    bw = zeros(length(time))

    return  StructArray{Comrade.EHTVisibilityDatum{Float64}}(
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
    obsamps = obs.amp::PyObject
    uamp = deepcopy(get(obsamps, Vector{Float64}, "u"))
    vamp = deepcopy(get(obsamps, Vector{Float64}, "v"))
    erramp = deepcopy(get(obsamps, Vector{Float64}, "sigma"))
    amps = deepcopy(get(obsamps, Vector{Float64}, "amp"))
    t1 = Symbol.(deepcopy(get(obsamps, Vector{String}, "t1")))
    t2 = Symbol.(deepcopy(get(obsamps, Vector{String}, "t2")))
    baseline = tuple.(t1, t2)
    time = deepcopy(get(obsamps, Vector{Float64}, "time"))
    freq = fill(obs.rf, length(time))

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
    obscp = obs.cphase::PyObject
    u1 = deepcopy((get(obscp, Vector{Float64}, "u1")))
    v1 = deepcopy((get(obscp, Vector{Float64}, "v1")))
    u2 = deepcopy((get(obscp, Vector{Float64}, "u2")))
    v2 = deepcopy((get(obscp, Vector{Float64}, "v2")))
    u3 = deepcopy((get(obscp, Vector{Float64}, "u3")))
    v3 = deepcopy((get(obscp, Vector{Float64}, "v3")))
    cp = deg2rad.(deepcopy((get(obscp, Vector{Float64}, "cphase"))))
    errcp = deg2rad.(deepcopy((get(obscp, Vector{Float64}, "sigmacp"))))

    t1 = Symbol.(deepcopy((get(obscp, Vector{String}, "t1"))))
    t2 = Symbol.(deepcopy((get(obscp, Vector{String}, "t2"))))
    t3 = Symbol.(deepcopy((get(obscp, Vector{String}, "t3"))))
    baseline = tuple.(t1, t2, t3)
    time = deepcopy(get(obscp, Vector{Float64}, "time"))
    freq = fill(obs.rf, length(time))
    bw = zeros(length(time))

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
    obslcamp = obs.logcamp::PyObject
    u1 = get(obslcamp, Vector{Float64}, "u1")
    v1 = get(obslcamp, Vector{Float64}, "v1")
    u2 = get(obslcamp, Vector{Float64}, "u2")
    v2 = get(obslcamp, Vector{Float64}, "v2")
    u3 = get(obslcamp, Vector{Float64}, "u3")
    v3 = get(obslcamp, Vector{Float64}, "v3")
    u4 = get(obslcamp, Vector{Float64}, "u4")
    v4 = get(obslcamp, Vector{Float64}, "v4")
    camp = ((get(obslcamp, Vector{Float64}, "camp")))
    errcamp = ((get(obslcamp, Vector{Float64}, "sigmaca")))

    t1 = Symbol.(((get(obslcamp, Vector{String}, "t1"))))
    t2 = Symbol.(((get(obslcamp, Vector{String}, "t2"))))
    t3 = Symbol.(((get(obslcamp, Vector{String}, "t3"))))
    t4 = Symbol.(((get(obslcamp, Vector{String}, "t4"))))
    baseline = tuple.(t1, t2, t3, t4)
    time = (get(obslcamp, Vector{Float64}, "time"))
    freq = fill(obs.rf, length(time))
    bw = zeros(length(time))

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
    u = get(obs.data, "u")
    v = get(obs.data, "v")

    # get visibilities
    c11 = get(obs.data, "rrvis")
    c12 = get(obs.data, "rlvis")
    c21 = get(obs.data, "lrvis")
    c22 = get(obs.data, "llvis")

    cohmat = StructArray{SMatrix{2,2,eltype(c11), 4}}((c11, c21, c12, c22))

    # get uncertainties
    e11 = get(obs.data, "rrsigma")
    e12 = get(obs.data, "rlsigma")
    e21 = get(obs.data, "lrsigma")
    e22 = get(obs.data, "llsigma")

    errmat = StructArray{SMatrix{2,2,eltype(e11), 4}}((e11, e21, e12, e22))


    # get timestamps and frequencies
    time = get(obs.data, "time")
    freq = fill(obs.rf, length(time))

    # get baseline info
    t1 = get(obs.data, Vector{Symbol}, "t1")
    t2 = get(obs.data, Vector{Symbol}, "t2")
    baseline = tuple.(t1, t2)

    # provide the polarization basis info
    single_polbasis = (CirBasis(), CirBasis())
    polbasis = fill(single_polbasis,length(u))

    # prepare output
    output = StructArray{EHTCoherencyDatum{eltype(u),
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


function getradec(obs)::Tuple{Float64, Float64}
    return (float(obs.ra), float(obs.dec))
end

function getmjd(obs)::Int
    return Int(obs.mjd)
end

function getsource(obs)::Symbol
    return Symbol(obs.source)
end


"""
    extract_amp(obs)
Extracts the visibility amplitudes from an ehtim observation object.

Any valid keyword arguments to `add_amp` in ehtim can be passed through extract_amp.

Returns an EHTObservation with visibility amplitude data
"""
function extract_amp(obsc; kwargs...)
    obs = obsc.copy()
    obs.add_scans()
    obs.reorder_tarr_snr()
    obs.add_amp(;kwargs...)
    data = getampfield(obs)
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf
    angles = get_fr_angles(obs)
    tarr = make_array_table(obsc)
    scans = Table((start=obs.scans[:,1], stop = obs.scans[:,2]))
    ac = _arrayconfig(data, angles, tarr, scans, bw)
    return Comrade.EHTObservation(data = data, mjd = mjd,
                   ra = ra, dec= dec,
                   config = ac,
                   bandwidth=bw,
                   source = source,
    )
end

function get_fr_angles(ehtobs)
    el1 = get(ehtobs.unpack(["el1"],ang_unit="rad"),"el1")
    el2 = get(ehtobs.unpack(["el2"],ang_unit="rad"),"el2")

    # read parallactic angles for each station
    par1 = get(ehtobs.unpack(["par_ang1"],ang_unit="rad"),"par_ang1")
    par2 = get(ehtobs.unpack(["par_ang2"],ang_unit="rad"),"par_ang2")
    return (el1, el2), (par1, par2)
end

function make_array_table(obs)
    return Table(
        sites = collect(Symbol.(get(obs.tarr, "site"))),
        X     = collect(get(obs.tarr, "x")),
        Y     = collect(get(obs.tarr, "y")),
        Z     = collect(get(obs.tarr, "z")),
        SEFD1 = collect(get(obs.tarr, "sefdr")),
        SEFD2 = collect(get(obs.tarr, "sefdl")),
        fr_parallactic = collect(get(obs.tarr, "fr_par")),
        fr_elevation = collect(get(obs.tarr, "fr_elev")),
        fr_offset = deg2rad.(collect(get(obs.tarr, "fr_off"))),
    )
end


"""
    extract_vis(obs)
Extracts the complex visibilities from an ehtim observation object

This grabs the raw `data` object from the obs object. Any keyword arguments are ignored.

Returns an EHTObservation with complex visibility data
"""
function extract_vis(obsc; kwargs...)
    obs = obsc.copy()
    obs.add_scans()
    obs.reorder_tarr_snr()

    data = getvisfield(obs)
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf
    angles = get_fr_angles(obs)
    tarr = make_array_table(obsc)
    scans = Table((start=obs.scans[:,1], stop = obs.scans[:,2]))
    ac = _arrayconfig(data, angles, tarr, scans, bw)
    return Comrade.EHTObservation(
                   data = data, mjd = mjd,
                   config=ac,
                   ra = ra, dec= dec,
                   bandwidth=bw,
                   source = source,
    )
end

"""
    extract_coherency(obs)
Extracts the coherency matrix from an ehtim observation object

This grabs the raw `data` object from the obs object. Any keyword arguments are ignored.

Returns an EHTObservation with coherency matrix
"""
function extract_coherency(obsc)
    obs = obsc.copy()
    obs.reorder_tarr_snr()
    obs.add_scans()
    data = getcoherency(obs)
    ra, dec = Comrade.getradec(obs)
    mjd = Comrade.getmjd(obs)
    source = Comrade.getsource(obs)
    bw = obs.bw
    rf = obs.rf
    angles = get_fr_angles(obs)
    tarr = make_array_table(obs)
    scans = Table((start=obs.scans[:,1], stop = obs.scans[:,2]))
    ac = _arrayconfig(data, angles, tarr, scans, bw)
    return Comrade.EHTObservation(data = data, mjd = mjd,
                   ra = ra, dec= dec,
                   config = ac,
                   bandwidth=bw,
                   source = source,
    )
end



function closure_designmat(scancp::Scan{A,B,C}, scanvis) where {A,B,C<:StructArray{<:EHTClosurePhaseDatum}}
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

function closure_designmat(scanlca::Scan{A,B,C}, scanvis) where {A,B,C<:StructArray{<:EHTLogClosureAmplitudeDatum}}
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
    dvis = extract_vis(obsc)
    st = scantable(dvis)

    minset, dmat = _minimal_closure(stlca, st)

    # now create EHTObservation
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf
    ac = arrayconfig(dvis)
    clac = ClosureConfig(dvis, dmat)
    return EHTObservation(data = minset, mjd = mjd,
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
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf


    cphase = EHTObservation(data = data, mjd = mjd,
                   config=nothing,
                   ra = ra, dec= dec,
                   bandwidth=bw,
                   source = source,
    )

    stcp = scantable(cphase)

    #Now make the vis obs
    dvis = extract_vis(obsc)
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

    ac = arrayconfig(dvis)
    clac = ClosureConfig(dvis, dmat)
    return  EHTObservation(data = data, mjd = mjd,
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
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf

    return EHTObservation(data = data, mjd = mjd,
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
    dvis = extract_vis(obsc)
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


    ac = arrayconfig(dvis)
    clac = ClosureConfig(dvis, dmat)
    return  EHTObservation(data = lcamp.data, mjd = lcamp.mjd,
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
    dvis = extract_vis(obsc)
    st = scantable(dvis)

    minset, dmat = _minimal_closure(stcp, st)
    # now create EHTObservation
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf
    ac = arrayconfig(dvis)
    clac = ClosureConfig(dvis, dmat)
    return EHTObservation(data = minset, mjd = mjd,
                          config=clac,
                          ra = ra, dec= dec,
                          bandwidth=bw,
                          source = source,
                        )
end

"""
    scan_average(obs; homogenize=true)

This homogenizes the scan times for an eht-imaging `Obsdata` object. This is needed
because eht-imaging has a bug that will sometimes create very small scans and
this can mess up both the closure construction and the gain scan times.
Note that this is only a problem if we
are fitting **scan averaged** data.
"""
function scan_average(obs)
    obsc = obs.copy()
    obsc.add_scans()
    obsc = obsc.avg_coherent(0.0, scan_avg=true)
    stimes = obsc.scans
    times = get(obsc.data, "time")
    @info "Before homogenizing we have $(length(unique(times))) unique times"

    for r in eachrow(stimes)
        sbegin, send = r
        indices = findall(x-> (sbegin < x <= send), times)
        times[indices] .= (send+sbegin)/2
    end
    set!(obsc.data, "time", times)
    @info "After homogenizing we have $(length(unique(times))) unique times"

    return obsc
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
function extract_cphase(obs; count="min-correct", cut_trivial=false, uvmin=0.1e9,  kwargs...)
    if count == "min-correct"
        return minimal_cphase(obs; cut_trivial, uvmin, kwargs...)
    else
        return _ehtim_cphase(obs; count=count, cut_trivial=cut_trivial, uvmin=uvmin, kwargs...)
    end
end



"""
    extract_lcamp(obs)
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
function extract_lcamp(obs; count="min-correct", kwargs...)
    if count == "min-correct"
        return minimal_lcamp(obs; kwargs...)
    else
        return _ehtim_lcamp(obs; count=count, kwargs...)
    end
end
