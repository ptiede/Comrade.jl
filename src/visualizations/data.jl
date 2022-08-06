export residuals, χ²

@recipe function f(m::AbstractModel, dvis::EHTObservation{T,A}; datamarker=:circle, datacolor=:grey) where {T,A<:EHTVisibilityDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "V (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :u)
    v = getdata(dvis, :v)
    uvdist = hypot.(u,v)
    vis = visibility.(dvis.data)
    error = getdata(dvis, :error)
    vre = real.(vis)
    vim = imag.(vis)
    #add data errorbars
    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := datacolor
        alpha := 0.5
        yerr := error
        linecolor := nothing
        label := "Data"
        uvdist, vre
    end

    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := :white
        markeralpha := 0.01
        markerstrokecolor := :black
        markerstrokealpha := 1.0
        linecolor :=nothing
        label := nothing
        yerr := error
        uvdist, vim
    end
    seriestype-->:scatter
    vmod = visibilities(m, u, v)
    labels --> "Model"
    uvdist, hcat(real.(vmod), imag.(vmod))
end

@recipe function f(dvis::EHTObservation{T,A};) where {T,A<:EHTVisibilityDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "V (Jy)"
    markershape --> :circle

    u = getdata(dvis, :u)
    v = getdata(dvis, :v)
    uvdist = hypot.(u,v)
    vis = visibility.(dvis.data)
    error = getdata(dvis, :error)
    vre = real.(vis)
    vim = imag.(vis)
    #add data errorbars
    @series begin
        seriestype := :scatter
        alpha := 0.5
        yerr := error
        linecolor := nothing
        label := "Real"
        uvdist, vre
    end

    @series begin
        seriestype := :scatter
        markeralpha := 0.1
        markerstrokecolor := :black
        markerstrokealpha := 1.0
        linecolor :=nothing
        label := nothing
        yerr := error
        label := "Imag"
        uvdist, vim
    end
end

@recipe function f(dvis::EHTObservation{T,A};) where {T,A<:EHTVisibilityAmplitudeDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "|V| (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :u)/1e9
    v = getdata(dvis, :v)/1e9
    uvdist = hypot.(u,v)
    amp = amplitude.(dvis.data)
    error = getdata(dvis, :error)
    #add data errorbars
    seriestype --> :scatter
    alpha --> 0.5
    yerr := error
    linecolor --> nothing
    label --> "Data"
    uvdist, amp
end

@recipe function f(acc::ArrayConfiguration)
    xguide --> "u (Gλ)"
    yguide --> "v (Gλ)"
    markershape --> :circle

    u, v = getuv(acc)
    #add data errorbars
    seriestype --> :scatter
    linecolor --> nothing
    aspect_ratio --> :equal
    label -->"Data"
    title --> "Frequency: $(acc.frequency/1e9) GHz"
    vcat(u/1e9,-u/1e9), vcat(v/1e9,-v/1e9)
end

@recipe function f(m::AbstractModel, dvis::EHTObservation{T,A}; datamarker=:circle, datacolor=:grey) where {T,A<:EHTVisibilityAmplitudeDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "V (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :u)
    v = getdata(dvis, :v)
    uvdist = hypot.(u,v)
    amp = amplitude.(dvis.data)
    error = getdata(dvis, :error)
    #add data errorbars
    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := datacolor
        alpha := 0.5
        yerr := error
        linecolor := nothing
        label := "Data"
        uvdist, amp
    end

    seriestype-->:scatter
    amod = amplitudes(m, u, v)
    labels --> "Model"
    uvdist, amod
end

export uvarea
function uvarea(d::EHTClosurePhaseDatum)
    u1, v1, u2, v2, u3, v3 = uvpositions(d)
    a = hypot(u1-u2, v1-v2)
    b = hypot(u2-u3, v2-v3)
    c = hypot(u3-u1, v3-v1)
    heron(a,b,c)
end

function heron(a,b,c)
    s = 0.5*(a+b+c)
    return sqrt(s*(s-a)*(s-b)*(s-c))
end

function uvarea(d::EHTLogClosureAmplitudeDatum)
    u1, v1, u2, v2, u3, v3, u4, v4 = uvpositions(d)
    a = hypot(u1-u2, v1-v2)
    b = hypot(u2-u3, v2-v3)
    c = hypot(u3-u4, v3-v4)
    d = hypot(u4-u1, v4-v1)
    h = hypot(u1-u3, v1-v3)
    return heron(a,b,h)+heron(c,d,h)
end

@recipe function f(dlca::EHTObservation{T,A}) where {T,A<:EHTLogClosureAmplitudeDatum}
    xguide --> "√(quadrangle area) (λ)"
    yguide --> "Log Clos. Amp."
    markershape --> :diamond
    area = sqrt.(uvarea.(dlca.data))
    phase = getdata(dlca, :amp)
    error = getdata(dlca, :error)
    #add data errorbars
    seriestype --> :scatter
    alpha := 0.5
    yerr := error
    label --> "Data"
    linecolor --> nothing
    area, phase
end


@recipe function f(m::AbstractModel, dlca::EHTObservation{T,A}; datamarker=:circle, datacolor=:grey) where {T,A<:EHTLogClosureAmplitudeDatum}
    xguide --> "√(quadrangle area) (λ)"
    yguide --> "Log Clos. Amp."
    markershape --> :diamond
    u1 = getdata(dlca, :u1)
    v1 = getdata(dlca, :v1)
    u2 = getdata(dlca, :u2)
    v2 = getdata(dlca, :v2)
    u3 = getdata(dlca, :u3)
    v3 = getdata(dlca, :v3)
    u4 = getdata(dlca, :u4)
    v4 = getdata(dlca, :v4)
    area = sqrt.(uvarea.(dlca.data))
    phase = getdata(dlca, :amp)
    error = getdata(dlca, :error)
    #add data errorbars
    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := datacolor
        alpha := 0.5
        yerr := error
        linecolor := nothing
        label := "Data"
        area, phase
    end

    seriestype-->:scatter
    amod = logclosure_amplitudes(m, u1, v1, u2, v2, u3, v3, u4, v4)
    labels --> "Model"
    area,amod
end

@recipe function f(dcp::EHTObservation{T,A}) where {T, A<:EHTClosurePhaseDatum}
    xguide --> "√(triangle area) (λ)"
    yguide --> "Phase (rad)"
    markershape --> :circle
    u1 = getdata(dcp, :u1)
    v1 = getdata(dcp, :v1)
    u2 = getdata(dcp, :u2)
    v2 = getdata(dcp, :v2)
    u3 = getdata(dcp, :u3)
    v3 = getdata(dcp, :v3)
    area = sqrt.(uvarea.(dcp.data))
    phase = getdata(dcp, :phase)
    error = getdata(dcp, :error)
    seriestype := :scatter
    alpha --> 0.5
    yerr := error
    linecolor --> nothing
    label --> "Data"
    area, phase
end

@recipe function f(m::AbstractModel, dcp::EHTObservation{T,A}; datamarker=:circle, datacolor=:grey) where {T,A<:EHTClosurePhaseDatum}
    xguide --> "√(triangle area) (λ)"
    yguide --> "Phase (rad)"
    markershape --> :diamond
    u1 = getdata(dcp, :u1)
    v1 = getdata(dcp, :v1)
    u2 = getdata(dcp, :u2)
    v2 = getdata(dcp, :v2)
    u3 = getdata(dcp, :u3)
    v3 = getdata(dcp, :v3)
    area = sqrt.(uvarea.(dcp.data))
    phase = getdata(dcp, :phase)
    error = getdata(dcp, :error)
    #add data errorbars
    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := datacolor
        alpha := 0.5
        yerr := error
        linecolor := nothing
        label --> "Data"
        area, phase
    end

    seriestype-->:scatter
    amod = closure_phases(m, u1, v1, u2, v2, u3, v3)
    labels --> "Model"
    area,amod
end

@userplot Residual

@recipe function f(h::Residual)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: AbstractModel) ||
        !(typeof(h.args[2]) <: EHTObservation)
        error("Residual should be given a model and data product.  Got: $(typeof(h.args))")
    end
    m, damp = h.args

    uvdist, res = residuals(m, damp)
    chi2 = sum(abs2, res)
    xguide --> "uv-distance (λ)"
    yguide --> "Normalized Residual"
    markershape --> :circle
    linecolor --> nothing
    legend --> false

    title --> @sprintf "<χ²> = %.2f" chi2/length(damp)
    uvdist, res
end


function χ²(m, data::EHTObservation)
    return sum(abs2, last(residuals(m, data)))
end

function χ²(m, data::EHTObservation...)
    return mapreduce(d->χ²(m, d), +, data)
end


function residuals(m, damp::EHTObservation{T, A}) where {T, A<:EHTVisibilityDatum}
    u = getdata(damp, :u)
    v = getdata(damp, :v)
    vis = StructArray{Complex{Float64}}((damp[:visr], damp[:visi]))
    mvis = visibilities(m, u, v)
    res = (vis - mvis)./getdata(damp, :error)
    re = real.(res)
    im = imag.(res)
    return hypot.(u, v), hcat(re, im)
end

function residuals(m, damp::EHTObservation{T, A}) where {T, A<:EHTVisibilityAmplitudeDatum}
    amp = getdata(damp, :amp)
    u = getdata(damp, :u)
    v = getdata(damp, :v)

    mamp = amplitudes(m, u, v)
    res = (amp - mamp)./getdata(damp, :error)
    return hypot.(u, v), res
end


function residuals(m, dcp::EHTObservation{T, A}) where {T, A<:EHTClosurePhaseDatum}
    area = sqrt.(uvarea.(dcp.data))
    phase = getdata(dcp, :phase)
    error = getdata(dcp, :error)

    mphase = closure_phases(m, dcp.config)
    res = zeros(length(phase))
    for i in eachindex(res)
        #s,c  = sincos(phase[i] - mphase[i])
        #dθ = atan(s,c)
        res[i] = abs(cis(phase[i]) - cis(mphase[i]))/error[i]
    end
return area, res
end


function residuals(m, dlca::EHTObservation{T, A}) where {T, A<:EHTLogClosureAmplitudeDatum}
    area = sqrt.(uvarea.(dlca.data))
    phase = getdata(dlca, :amp)
    error = getdata(dlca, :error)

    mphase = logclosure_amplitudes(m, dlca.config)
    res = (phase- mphase)./error
    return area, res
end
