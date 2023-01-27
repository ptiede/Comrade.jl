export extract_amp, extract_vis, extract_lcamp, extract_cphase,
       extract_coherency,
       load_ehtim_uvfits, scan_average
using PyCall: set!, PyObject


struct EHTIMObs{O} <: Observation{Float64}
    obsdata::O
end

function getproperty(obs::EHTIMObs, v::Any)
    out = getproperty(getfield(obs, :obsdata), v)
    if out isa PyObject
        return EHTIMObs(out)
    else
        return out
    end
end
obsdata(obs::EHTIMObs) = getfield(obs, :obsdata)
scans(obs::EHTIMObs)         = Table(start=@view obsdata(obs).scans[:,1], stop=@view obsdata(obs).scans[:,2])
telescope_array(d::EHTIMObs) = make_array_table(obsdata(d))
bandwidth(d::EHTIMObs)       = obsdata(bw)
source(d::EHTIMObs)          = Symbol(obsdata(d).source)
ra(d::EHTIMObs)              = obsdata(d).ra
dec(d::EHTIMObs)             = obsdata(d).dec



ehtimstokes(::IPol) = ""
ehtimstokes(::QPol) = "q"
ehtimstokes(::UPol) = "u"
ehtimstokes(::VPol) = "v"

function getvisfield(obs::EHTIMObs, pol::StokesBasis = IPol())
    config = ehtarrayconfig(obs, pol)
    e = ehtimstokes(pol)
    err = deepcopy((get(obsamps, Vector{Float64}, e*"sigma")))
    vis = deepcopy((get(obsamps, Vector{Complex{Float64}}, e*"vis")))
    return  config, vis, err
end


function getampfield(obs, pol::StokesBasis = IPol(), debias=false)

    return config, amp, err
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
function extract_amp(obs::EHTIMObs, pol=IPol(); debias=false, kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()

    config = ehtarrayconfig(obs, pol)

    e = ehtimstokes(pol)
    err = deepcopy((get(obsamps, Vector{Float64}, e*"sigma")))
    vis = deepcopy((get(obsamps, Vector{Complex{Float64}}, e*"vis")))

    if debias
        amp = sqrt.(max.(abs.(vis).^2 - err.^2, 0.0))
    else
        amp = abs.(vis)
    end

    return EHTDataTable{EHTAmplitudeDatum}(
                 meas, error, config,
                 mjd(obs), ra(obs), dec(obs),
                 source(obs), obs
                )
end

function get_fr_angles(obs::EHTIMObs)
    ehtobs = obsdata(obs)
    el1 = get(ehtobs.unpack(["el1"],ang_unit="rad"),"el1")
    el2 = get(ehtobs.unpack(["el2"],ang_unit="rad"),"el2")

    # read parallactic angles for each station
    par1 = get(ehtobs.unpack(["par_ang1"],ang_unit="rad"),"par_ang1")
    par2 = get(ehtobs.unpack(["par_ang2"],ang_unit="rad"),"par_ang2")
    return (el1, el2), (par1, par2)
end

function make_array_table(obse::EHTIMObs)
    obs = obsdata(obse)
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
function extract_vis(obsc::EHTIMObs, pol=IPol(); kwargs...)
    obs = obsc.copy()
    obs.reorder_tarr_snr()

    config = ehtarrayconfig(obs, pol)
    e = ehtimstokes(pol)
    error = deepcopy((get(obsdata(obs), Vector{Float64}, e*"sigma")))
    meas = deepcopy((get(obsdata(obs), Vector{Complex{Float64}}, e*"vis")))

    return EHTDataTable{EHTVisibilityDatum}(
                meas, error, config,
                mjd(obs), ra(obs), dec(obs),
                source(obs), obs
                )
end

"""
    extract_coherency(obs)
Extracts the coherency matrix from an ehtim observation object

This grabs the raw `data` object from the obs object. Any keyword arguments are ignored.

Returns an EHTObservation with coherency matrix
"""
function extract_coherency(obsc::EHTIMObs)
    obs = obsc.copy()
    obs.reorder_tarr_snr()
    obsd = obsdata(obs)

    config = ehtarrayconfig(obs, pol)
    # get visibilities
    c11 = get(obsd.data, "rrvis")
    c12 = get(obsd.data, "rlvis")
    c21 = get(obsd.data, "lrvis")
    c22 = get(obsd.data, "llvis")

    meas = StructArray{SMatrix{2,2,eltype(c11), 4}}((c11, c21, c12, c22))

    # get uncertainties
    e11 = get(obsd.data, "rrsigma")
    e12 = get(obsd.data, "rlsigma")
    e21 = get(obsd.data, "lrsigma")
    e22 = get(obsd.data, "llsigma")

    error = StructArray{SMatrix{2,2,eltype(e11), 4}}((e11, e21, e12, e22))

    return EHTDataTable{EHTCoherencyDatum}(
                meas, error, config,
                mjd(obs), ra(obs), dec(obs),
                source(obs), obs
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


function _ehtim_cphase(obs; count="max", cut_trivial=false, uvmin=0.1e9, kwargs...)
    obs = obsc.copy()

    # cut 0 baselines since these are trivial triangles

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
    clac = ClosureConfig(ac, dmat)
    return  EHTObservation(data = data, mjd = mjd,
                           config=clac,
                           ra = ra, dec= dec,
                           bandwidth=bw,
                           source = source
                          )
end

function _make_lcamp(obsc, count="max"; kwargs...)
    obs = obsc.copy()

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
function scan_average(obs::EHTIMObs)
    obsc = copy(obs.obsdata)
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

function add_fractional_noise(obs::EHTIMObs, frac_noise::Real, debias::Bool = false)
    obsdatac = obs.obsdata.add_fractional(frac_noise, debias)
    return EHTIMObs(obsdatac)
end

function coherency_average(obs::EHTIMObs, dt::Real)
    obsavg = obs.obsdata.avg_coherent(dt)
    return EHTIMObs(obsavg)
end

function extract_table(obs::EHTIMObs, ::ComplexVis; kwargs...)
    return extract_vis(obs.obsdata, kwargs...)
end

function extract_table(obs::EHTIMObs, p::CPhase; kwargs...)
    if p.minimal
        count = "min-correct"
    else
        count = "max"
    end
    dcp = extract_cphase(obs.obsdata, count=count, kwargs...)
end

function extract_table(obs::EHTIMObs, p::LogCamp; kwargs...)
    if p.minimal
        count = "min-correct"
    else
        count = "max"
    end
    dlca = extract_cphase(obs.obsdata, count=count, kwargs...)
end


function extract_table(obs::EHTIMObs, ::Coherency; kwargs...)
    return extract_coherency(obs.obsdata)
end

function extract_table(obs::EHTIMObs, ::Amplitude, kwargs...)
    return extract_amp(obs.obsdata, kwargs...)
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
