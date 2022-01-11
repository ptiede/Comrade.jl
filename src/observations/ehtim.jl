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

"""
    extract_cphase(obs)
Extracts the closure phases from an ehtim observation object

Returns an EHTObservation with closure phases datums
"""
function extract_cphase(obs; kwargs...)
    obs.add_cphase(;kwargs...)
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

end


"""
    extract_lcamp(obs)
Extracts the log-closure amp. from an ehtim observation object

Returns an EHTObservation with closure amp. datums
"""
function extract_lcamp(obs; kwargs...)
    obs.add_logcamp(;kwargs...)
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
