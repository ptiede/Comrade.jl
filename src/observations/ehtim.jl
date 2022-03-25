export extract_amp, extract_vis, extract_lcamp, extract_cphase

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
    freq = zeros(length(time))
    bw = zeros(length(time))

    return  StructArray{Comrade.EHTVisibilityDatum{Float64}}(
        visr = real.(vis),
        visi = imag.(vis),
        u = u,
        v = v,
        error = err,
        time = time,
        frequency = freq,
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
    freq = zeros(length(time))
    bw = zeros(length(time))

    return  StructArray{Comrade.EHTVisibilityAmplitudeDatum{Float64}}(
        amp = amps,
        u = uamp,
        v = vamp,
        error = erramp,
        time = time,
        frequency = freq,
        bandwidth = bw,
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
    freq = zeros(length(time))
    bw = zeros(length(time))

    return StructArray{Comrade.EHTClosurePhaseDatum{Float64}}(
        phase = cp,
        u1 = u1,
        v1 = v1,
        u2 = u2,
        v2 = v2,
        u3 = u3,
        v3 = v3,
        error = errcp,
        time = time,
        frequency = freq,
        bandwidth = bw,
        triangle = baseline
    )

end

function getlcampfield(obs)
    obslcamp = obs.logcamp::PyObject
    u1 = deepcopy((get(obslcamp, Vector{Float64}, "u1")))
    v1 = deepcopy((get(obslcamp, Vector{Float64}, "v1")))
    u2 = deepcopy((get(obslcamp, Vector{Float64}, "u2")))
    v2 = deepcopy((get(obslcamp, Vector{Float64}, "v2")))
    u3 = deepcopy((get(obslcamp, Vector{Float64}, "u3")))
    v3 = deepcopy((get(obslcamp, Vector{Float64}, "v3")))
    u4 = deepcopy((get(obslcamp, Vector{Float64}, "u4")))
    v4 = deepcopy((get(obslcamp, Vector{Float64}, "v4")))
    camp = deepcopy((get(obslcamp, Vector{Float64}, "camp")))
    errcamp = deepcopy((get(obslcamp, Vector{Float64}, "sigmaca")))

    t1 = Symbol.(deepcopy((get(obslcamp, Vector{String}, "t1"))))
    t2 = Symbol.(deepcopy((get(obslcamp, Vector{String}, "t2"))))
    t3 = Symbol.(deepcopy((get(obslcamp, Vector{String}, "t3"))))
    t4 = Symbol.(deepcopy((get(obslcamp, Vector{String}, "t4"))))
    baseline = tuple.(t1, t2, t3, t4)
    time = deepcopy(get(obslcamp, Vector{Float64}, "time"))
    freq = zeros(length(time))
    bw = zeros(length(time))

    return StructArray{Comrade.EHTLogClosureAmplitudeDatum{Float64}}(
        amp = camp,
        u1 = u1,
        v1 = v1,
        u2 = u2,
        v2 = v2,
        u3 = u3,
        v3 = v3,
        u4 = u4,
        v4 = v4,
        error = errcamp,
        time = time,
        frequency = freq,
        bandwidth = bw,
        quadrangle = baseline
    )

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
Extracts the visibility amplitudes from an ehtim observation object

Returns an EHTObservation with visibility amplitude data
"""
function extract_amp(obs; kwargs...)
    obs.add_amp(;kwargs...)
    data = getampfield(obs)
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf

    return Comrade.EHTObservation(data = data, mjd = mjd,
                   ra = ra, dec= dec,
                   bandwidth=bw, frequency=rf,
                   source = source,
    )
end


"""
    extract_vis(obs)
Extracts the complex visibilities from an ehtim observation object

Returns an EHTObservation with complex visibility data
"""
function extract_vis(obs)
    data = getvisfield(obs)
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf
    return Comrade.EHTObservation(
                   data = data, mjd = mjd,
                   ra = ra, dec= dec,
                   bandwidth=bw, frequency=rf,
                   source = source,
    )
end


function closure_designmat(scancp::Scan{A,B,C}, scanvis) where {A,B,C<:StructArray{<:EHTClosurePhaseDatum}}
    ant1, ant2, ant3 = baselines(scancp)
    antvis1, antvis2 = baselines(scanvis)
    design_mat = zeros(length(scancp), length(scanvis))
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
        (((av1 == a1)&(av2 == a2)) | (av1 == a2)&(av2 == a1))&& (design_mat[j,i] = 1.0)

        #check leg 2
        (((av1 == a3)&(av2 == a4)) | (av1 == a4)&(av2 == a3))&& (design_mat[j,i] = 1.0)

        #check leg 3
        (((av1 == a1)&(av2 == a4)) | (av1 == a4)&(av2 == a1))&& (design_mat[j,i] = -1.0)

        #check leg 4
        (((av1 == a2)&(av2 == a3)) | (av1 == a3)&(av2 == a2))&& (design_mat[j,i] = -1.0)

    end
    return design_mat
end

function minimal_lcamp(obsc; kwargs...)

    obs = obsc.copy()
    # reorder to maximize the snr
    obs.reorder_tarr_snr()

    lcamp = extract_lcamp(obs; count="max")
    stlca = scantable(lcamp)

    #Now make the vis obs
    dvis = extract_vis(obs)
    st = scantable(dvis)

    minset = _minimal_closure(stlca, st)

    # now create EHTObservation
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf
    return EHTObservation(data = minset, mjd = mjd,
                          ra = ra, dec= dec,
                          bandwidth=bw, frequency=rf,
                          source = source
                        )

end

function minimal_cphase(obsc; kwargs...)
    # compute a maximum set of closure phases
    obs = obsc.copy()

    # reorder to maximize the snr
    obs.reorder_tarr_snr()


    cphase = extract_cphase(obs; count="max")
    stcp = scantable(cphase)

    #Now make the vis obs
    dvis = extract_vis(obs)
    st = scantable(dvis)

    minset = _minimal_closure(stcp, st)
    # now create EHTObservation
    ra, dec = getradec(obs)
    mjd = getmjd(obs)
    source = getsource(obs)
    bw = obs.bw
    rf = obs.rf
    return EHTObservation(data = minset, mjd = mjd,
                          ra = ra, dec= dec,
                          bandwidth=bw, frequency=rf,
                          source = source
                        )
end

function _minimal_closure(stcl, st)

    # Determine the number of timestamps containing a closure triangle
    N_times_cl = length(stcl)
    T = typeof(stcl[1][1])
    # loop over all timestamps
    minset = StructVector{T}(undef, 0)
    for i in 1:N_times_cl

        # Get our clhase scan
        scancl = stcl[i]

        # sort by clhase SNR to we form a nice minimal set
        snr = inv.(scancl.scan.error)
        ind_snr = sortperm(snr)
        scancl = scancl[ind_snr]
        snr = snr[ind_snr]

        # now find the visibility scan that matches the closure one based on timestamps
        ind_here_bl = findfirst(==(scancl.time), st.times)
        scanvis = st[ind_here_bl]


        # initialize the design matrix
        design_mat = closure_designmat(scancl, scanvis)

        # determine the expected size of the minimal set
        # this is needed to make sure we aren't killing too many triangles
        nmin = rank(design_mat)

        # print some info
        #println("For timestamp $(stcl.times[i]):")

        # get the current stations
        #println("Observing stations are $(stations(scancl))")

        println("Size of maximal set of closure products = $(length(scancl))")
        println("Size of minimal set of closure products = $(nmin)")
        println("...")

        ##########################################################
        # start of loop to recover minimal set
        ##########################################################
        minset_scan = minimal_closure_scan(scancl, scanvis, nmin)
        minset = append!!(minset, minset_scan.scan)
    end
    return minset
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
    if length(scankeep) != nmin
        @error "minimal set not found $(length(scankeep)) $(nmin)"
        throw("No minimal set found at time $(scancl.time)")
    end
    return scankeep
end

"""
    extract_cphase(obs)
Extracts the closure phases from an ehtim observation object

Returns an EHTObservation with closure phases datums
"""
function extract_cphase(obs; count="min", kwargs...)

    if count == "min"
        return minimal_cphase(obs; kwargs...)
    else
        obs.add_cphase(;count, kwargs...)
        data = getcpfield(obs)
        ra, dec = getradec(obs)
        mjd = getmjd(obs)
        source = getsource(obs)
        bw = obs.bw
        rf = obs.rf

        return Comrade.EHTObservation(data = data, mjd = mjd,
                       ra = ra, dec= dec,
                       bandwidth=bw, frequency=rf,
                       source = source,
        )
    #end

end


"""
    extract_lcamp(obs)
Extracts the log-closure amp. from an ehtim observation object

Returns an EHTObservation with closure amp. datums
"""
function extract_lcamp(obs; count="min", kwargs...)
    if count == "min"
        return minimal_cphase(obs; kwargs...)
    else

        obs.add_logcamp(;count, kwargs...)
        data = getlcampfield(obs)
        ra, dec = getradec(obs)
        mjd = getmjd(obs)
        source = getsource(obs)
        bw = obs.bw
        rf = obs.rf

        return Comrade.EHTObservation(data = data, mjd = mjd,
                       ra = ra, dec= dec,
                       bandwidth=bw, frequency=rf,
                       source = source,
        )
    end
end
