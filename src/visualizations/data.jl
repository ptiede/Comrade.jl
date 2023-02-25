export residuals, chi2

@recipe function f(m::AbstractModel, dvis::EHTObservation{T,A}; datamarker=:circle, datacolor=:grey) where {T,A<:EHTVisibilityDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "V (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
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
    vmod = visibilitymap(m, arrayconfig(dvis))
    labels --> "Model"
    uvdist, hcat(real.(vmod), imag.(vmod))
end

@recipe function f(m::AbstractModel, dvis::EHTObservation{T,A}; datamarker=:circle, datacolor=:grey) where {T,A<:EHTCoherencyDatum}

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    uvdist = hypot.(u,v)./1e9
    vis = dvis[:measurement]
    error = getdata(dvis, :error)
    vmod = visibilitymap(m, (U=u, V=v, T=dvis[:T], F=dvis[:F]))
    layout := (2,2)
    #add data errorbars
    @series begin
        v = getindex.(vis, 1, 1)
        vm = getindex.(vmod, 1, 1)
        err = getindex.(error, 1, 1)
        subplot := 1
        yguide := "C₁₁ (Jy)"
        seriestype := :scatter

        @series begin
            markershape := [:+ :x]
            markercolor := datacolor
            seriestype := :scatter
            alpha := 0.9
            yerr := err
            linecolor := nothing
            label := ["Data Real" "DataImag"]
            uvdist, hcat(real.(v), imag.(v))
        end

        label := ["Model Real" "Model Imag"]
        markershape := :o
        markercolor := [:red :orange]
        markerstrokecolor := [:red :orange]
        markeralpha := 0.2
        uvdist, hcat(real.(vm), imag.(vm))
    end

    @series begin
        subplot := 2
        v = getindex.(vis, 1, 2)
        vm = getindex.(vis, 1, 2)
        err = getindex.(error, 1, 2)
        yguide := "C₁₂ (Jy)"
        legend := nothing
        seriestype := :scatter
        @series begin
            markershape := [:+ :x]
            markercolor := datacolor
            seriestype := :scatter
            alpha := 0.9
            yerr := err
            linecolor := nothing
            label := ["Data Real" "DataImag"]
            uvdist, hcat(real.(v), imag.(v))
        end

        label := ["Model Real" "Model Imag"]
        markershape := :o
        markercolor := [:red :orange]
        markerstrokecolor := [:red :orange]
        markeralpha := 0.2
        uvdist, hcat(real.(vm), imag.(vm))
    end


    @series begin
        subplot := 3
        v = getindex.(vis, 2, 1)
        vm = getindex.(vis, 2, 1)
        err = getindex.(error, 2, 1)
        legend := nothing
        seriestype := :scatter
        yguide := "C₂₁ (Jy)"
        @series begin
            markershape := [:+ :x]
            markercolor := datacolor
            seriestype := :scatter
            alpha := 0.9
            yerr := err
            linecolor := nothing
            label := ["Data Real" "DataImag"]
            uvdist, hcat(real.(v), imag.(v))
        end

        xguide := "uv distance (Gλ)"
        label := ["Model Real" "Model Imag"]
        markershape := :o
        markercolor := [:red :orange]
        markerstrokecolor := [:red :orange]
        markeralpha := 0.2
        uvdist, hcat(real.(vm), imag.(vm))
    end

    @series begin
        subplot := 4
        v = getindex.(vis, 2, 2)
        vm = getindex.(vis, 2, 2)
        err = getindex.(error, 2, 2)
        legend := nothing
        seriestype := :scatter
        yguide := "C₂₂ (Jy)"
        @series begin
            markershape := [:+ :x]
            markercolor := datacolor
            seriestype := :scatter
            alpha := 0.9
            yerr := err
            linecolor := nothing
            label := ["Data Real" "DataImag"]
            uvdist, hcat(real.(v), imag.(v))
        end

        label := ["Model Real" "Model Imag"]
        markershape := :o
        markercolor := [:red :orange]
        markerstrokecolor := [:red :orange]
        markeralpha := 0.2
        uvdist, hcat(real.(vm), imag.(vm))

    end

end

@recipe function f(dvis::EHTObservation{T,A};) where {T,A<:EHTVisibilityDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "V (Jy)"
    markershape --> :circle

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
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

@recipe function f(dvis::EHTObservation{T,A};) where {T,A<:EHTCoherencyDatum}
    markershape --> :circle

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    uvdist = hypot.(u,v)
    coh = dvis[:measurement]
    error = dvis[:error]
    layout := (2,2)

    #add data errorbars
    @series begin
        xguide --> "uv-distance (λ)"
        subplot := 1
        yguide := "RR (Jy)"
        vre = real.(getindex.(coh, 1, 1))
        vim = imag.(getindex.(coh, 1, 1))
        err = getindex.(error, 1, 1)
        @series begin
            seriestype := :scatter
            alpha := 0.5
            yerr := err
            linecolor := nothing
            label := "Real"
            uvdist, vre
        end

            seriestype := :scatter
            markeralpha := 0.5
            markerstrokecolor := :black
            markerstrokealpha := 1.0
            linecolor :=nothing
            label := nothing
            yerr := err
            label := "Imag"
            uvdist, vim
    end

    @series begin
        xguide --> "uv-distance (λ)"
        subplot := 2
        yguide := "RL (Jy)"
        vre = real.(getindex.(coh, 1, 2))
        vim = imag.(getindex.(coh, 1, 2))
        err = getindex.(error, 1, 2)
        @series begin
            seriestype := :scatter
            alpha := 0.5
            yerr := err
            linecolor := nothing
            label := "Real"
            uvdist, vre
        end

            seriestype := :scatter
            markeralpha := 0.5
            markerstrokecolor := :black
            markerstrokealpha := 1.0
            linecolor :=nothing
            legend := nothing
            label := nothing
            yerr := err
            label := "Imag"
            uvdist, vim
    end

    @series begin
        xguide --> "uv-distance (λ)"
        subplot := 3
        yguide := "LR (Jy)"
        vre = real.(getindex.(coh, 2, 1))
        vim = imag.(getindex.(coh, 2, 1))
        err = getindex.(error, 1, 1)
        @series begin
            seriestype := :scatter
            alpha := 0.5
            yerr := err
            linecolor := nothing
            label := "Real"
            uvdist, vre
        end

            seriestype := :scatter
            markeralpha := 0.5
            markerstrokecolor := :black
            markerstrokealpha := 1.0
            linecolor :=nothing
            label := nothing
            legend := nothing
            yerr := err
            label := "Imag"
            uvdist, vim
    end
    @series begin
        xguide --> "uv-distance (λ)"
        subplot := 4
        yguide := "LL (Jy)"
        vre = real.(getindex.(coh, 2, 2))
        vim = imag.(getindex.(coh, 2, 2))
        err = getindex.(error, 1, 1)
        @series begin
            seriestype := :scatter
            alpha := 0.5
            yerr := err
            linecolor := nothing
            label := "Real"
            uvdist, vre
        end

            seriestype := :scatter
            markeralpha := 0.5
            markerstrokecolor := :black
            markerstrokealpha := 1.0
            linecolor :=nothing
            label := nothing
            legend := nothing
            yerr := err
            label := "Imag"
            uvdist, vim
    end

end

@recipe function f(dvis::EHTObservation{T,A};) where {T,A<:EHTVisibilityAmplitudeDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "|V| (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :U)/1e9
    v = getdata(dvis, :V)/1e9
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
    title --> "Frequency: $(first(acc.data.F)/1e9) GHz"
    vcat(u/1e9,-u/1e9), vcat(v/1e9,-v/1e9)
end

@recipe function f(m::AbstractModel, dvis::EHTObservation{T,A}; datamarker=:circle, datacolor=:grey) where {T,A<:EHTVisibilityAmplitudeDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "V (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
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
    amod = abs.(visibilitymap(m, arrayconfig(dvis)))
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
    phase = getdata(dlca, :measurement)
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
    u1 = getdata(dlca, :U1)
    v1 = getdata(dlca, :V1)
    u2 = getdata(dlca, :U2)
    v2 = getdata(dlca, :V2)
    u3 = getdata(dlca, :U3)
    v3 = getdata(dlca, :V3)
    u4 = getdata(dlca, :U4)
    v4 = getdata(dlca, :V4)
    area = sqrt.(uvarea.(dlca.data))
    phase = getdata(dlca, :measurement)
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
    amod = logclosure_amplitudes(m, dlca.config)
    labels --> "Model"
    area,amod
end

@recipe function f(dcp::EHTObservation{T,A}) where {T, A<:EHTClosurePhaseDatum}
    xguide --> "√(triangle area) (λ)"
    yguide --> "Phase (rad)"
    markershape --> :circle
    u1 = getdata(dcp, :U1)
    v1 = getdata(dcp, :V1)
    u2 = getdata(dcp, :U2)
    v2 = getdata(dcp, :V2)
    u3 = getdata(dcp, :U3)
    v3 = getdata(dcp, :V3)
    area = sqrt.(uvarea.(dcp.data))
    phase = getdata(dcp, :measurement)
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
    u1 = getdata(dcp, :U1)
    v1 = getdata(dcp, :V1)
    u2 = getdata(dcp, :U2)
    v2 = getdata(dcp, :V2)
    u3 = getdata(dcp, :U3)
    v3 = getdata(dcp, :V3)
    area = sqrt.(uvarea.(dcp.data))
    phase = getdata(dcp, :measurement)
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
    amod = closure_phases(m, dcp.config)
    labels --> "Model"
    area,amod
end

@userplot Residual

export ndata
ndata(d::EHTObservation) = length(d)
ndata(d::EHTObservation{T, D}) where {T, D<:EHTVisibilityDatum} = 2*length(d)
ndata(d::EHTObservation{T, D}) where {T, D<:EHTCoherencyDatum} = 8*length(d)

@recipe function f(h::Residual)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: AbstractModel) ||
        !(typeof(h.args[2]) <: EHTObservation)
        error("Residual should be given a model and data product.  Got: $(typeof(h.args))")
    end
    m, damp = h.args
    uvdist, res = residuals(m, damp)
    c2 = chi2(m, damp)
    # title-->"Norm. Residuals"
    legend-->nothing

    if damp isa EHTObservation{<:Any, <:EHTCoherencyDatum}
        layout := (2,2)
        res2 = reinterpret(reshape, Float64, res)'
        @series begin
            yguide := "RR"
            subplot := 1
            seriestype := :scatter
            alpha := 0.5
            linecolor := nothing
            title --> @sprintf "χ² = %.2f" sum(abs2, filter(!isnan, @view res2[:,1:2]))/(2*size(res2,1))
            uvdist, res2[:,1:2]
        end
        @series begin
            xguide --> "uv-distance (λ)"
            yguide := "LR"
            subplot := 3
            seriestype := :scatter
            alpha := 0.5
            linecolor := nothing
            title --> @sprintf "χ² = %.2f" sum(abs2, filter(!isnan, @view res2[:,3:4]))/(2*size(res2,1))
            uvdist, res2[:,3:4]
        end
        @series begin
            yguide := "RL"
            subplot := 2
            seriestype := :scatter
            alpha := 0.5
            linecolor := nothing
            title --> @sprintf "χ² = %.2f" sum(abs2, filter(!isnan, @view res2[:,5:6]))/(2*size(res2,1))
            uvdist, res2[:,5:6]
        end
        @series begin
            yguide := "LL"
            subplot := 4
            seriestype := :scatter
            alpha := 0.5
            linecolor := nothing
            title --> @sprintf "χ² = %.2f" sum(abs2, filter(!isnan, @view res2[:,7:8]))/(2*size(res2,1))
            uvdist, res2[:,7:8]
        end
    else
        xguide --> "uv-distance (λ)"
        yguide --> "Normalized Residual"
        markershape --> :circle
        linecolor --> nothing
        legend --> false

        title --> @sprintf "<χ²> = %.2f" c2/ndata(damp)
        return uvdist, res
    end
end


function chi2(m, data::EHTObservation)
    return sum(x->abs2.(x), last(residuals(m, data)))
end

function chi2(m, data::EHTObservation{T, A}) where {T, A<:EHTCoherencyDatum}
    res = last(residuals(m, data))
    return mapreduce(+, 1:4) do i
        sum(abs2, filter(!isnan, getproperty(res, i)))
    end
end


function chi2(m, data::EHTObservation...)
    return mapreduce(d->chi2(m, d), +, data)
end


function residuals(m, dvis::EHTObservation{T, A}) where {T, A<:EHTVisibilityDatum}
    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    vis = dvis[:measurement]
    mvis = visibilitymap(m, (U=u, V=v, T=dvis[:T], F=dvis[:F]))
    res = (vis - mvis)./getdata(dvis, :error)
    re = real.(res)
    im = imag.(res)
    return hypot.(u, v), hcat(re, im)
end

function residuals(m, dvis::EHTObservation{T, A}) where {T, A<:EHTCoherencyDatum}
    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    coh = dvis[:measurement]
    mcoh = visibilitymap(m, (U=u, V=v, T=dvis[:T], F=dvis[:F]))
    res = map((x,y,z)->((x .- y)./z), coh, mcoh, dvis[:error])
    return hypot.(u, v), StructArray(res)
end


function residuals(m, damp::EHTObservation{T, A}) where {T, A<:EHTVisibilityAmplitudeDatum}
    amp = getdata(damp, :measurement)
    u = getdata(damp, :U)
    v = getdata(damp, :V)

    mamp = amplitudes(m, (U=u, V=v))
    res = (amp - mamp)./getdata(damp, :error)
    return hypot.(u, v), res
end


function residuals(m, dcp::EHTObservation{T, A}) where {T, A<:EHTClosurePhaseDatum}
    area = sqrt.(uvarea.(dcp.data))
    phase = getdata(dcp, :measurement)
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
    phase = getdata(dlca, :measurement)
    error = getdata(dlca, :error)

    mphase = logclosure_amplitudes(m, dlca.config)
    res = (phase- mphase)./error
    return area, res
end
