
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
        label := "Data"
        area, phase
    end

    seriestype-->:scatter
    amod = closure_phases(m, u1, v1, u2, v2, u3, v3)
    labels --> "Model"
    area,amod
end
