export residuals, chi2

@recipe function f(m::AbstractModel, dvis::EHTObservationTable{A}; datamarker=:circle, datacolor=:grey) where {A<:EHTVisibilityDatum}
    xguide --> "uv-distance (λ)"
    yguide --> "V (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    uvdist = hypot.(u,v)
    vis = dvis[:measurement]
    noise = getdata(dvis, :noise)
    vre = real.(vis)
    vim = imag.(vis)
    #add data noisebars
    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := datacolor
        alpha := 0.5
        yerr := noise
        linecolor := nothing
        label := "Data"
        uvdist./1e9, vre
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
        yerr := noise
        uvdist./1e9, vim
    end


    seriestype-->:scatter
    vmod = visibilitymap(m, arrayconfig(dvis))
    labels --> "Model"
    uvdist./1e9, hcat(real.(vmod), imag.(vmod))
end

@recipe function f(m::AbstractModel, dvis::EHTObservationTable{A}; datamarker=:circle, datacolor=:grey) where {A<:EHTCoherencyDatum}

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    uvdist = hypot.(u,v)./1e9
    vis = dvis[:measurement]
    noise = getdata(dvis, :noise)
    vmod = visibilitymap(m, (U=u, V=v, T=dvis[:T], F=dvis[:F]))
    layout := (2,2)
    #add data noisebars
    @series begin
        v = getindex.(vis, 1, 1)
        vm = getindex.(vmod, 1, 1)
        err = getindex.(noise, 1, 1)
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
        err = getindex.(noise, 1, 2)
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
        err = getindex.(noise, 2, 1)
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
        err = getindex.(noise, 2, 2)
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

@recipe function f(dvis::EHTObservationTable{A};) where {A<:EHTVisibilityDatum}
    xguide --> "uv-distance (Gλ)"
    yguide --> "V (Jy)"
    markershape --> :circle

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    uvdist = hypot.(u,v)
    vis = dvis[:measurement]
    noise = getdata(dvis, :noise)
    vre = real.(vis)
    vim = imag.(vis)
    #add data noisebars
    @series begin
        seriestype := :scatter
        alpha := 0.5
        yerr := noise
        linecolor := nothing
        label := "Real"
        uvdist./1e9, vre
    end

    @series begin
        seriestype := :scatter
        markeralpha := 0.1
        markerstrokecolor := :black
        markerstrokealpha := 1.0
        linecolor :=nothing
        label := nothing
        yerr := noise
        label := "Imag"
        uvdist./1e9, vim
    end
end

@recipe function f(dvis::EHTObservationTable{A};) where {A<:EHTCoherencyDatum}
    markershape --> :circle

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    uvdist = hypot.(u,v)
    coh = dvis[:measurement]
    noise = dvis[:noise]
    layout := (2,2)

    #add data noisebars
    @series begin
        xguide --> "uv-distance (Gλ)"
        subplot := 1
        yguide := "RR (Jy)"
        vre = real.(getindex.(coh, 1, 1))
        vim = imag.(getindex.(coh, 1, 1))
        err = getindex.(noise, 1, 1)
        @series begin
            seriestype := :scatter
            alpha := 0.5
            yerr := err
            linecolor := nothing
            label := "Real"
            uvdist./1e9, vre
        end

            seriestype := :scatter
            markeralpha := 0.5
            markerstrokecolor := :black
            markerstrokealpha := 1.0
            linecolor :=nothing
            label := nothing
            yerr := err
            label := "Imag"
            uvdist./1e9, vim
    end

    @series begin
        xguide --> "uv-distance (Gλ)"
        subplot := 2
        yguide := "RL (Jy)"
        vre = real.(getindex.(coh, 1, 2))
        vim = imag.(getindex.(coh, 1, 2))
        err = getindex.(noise, 1, 2)
        @series begin
            seriestype := :scatter
            alpha := 0.5
            yerr := err
            linecolor := nothing
            label := "Real"
            uvdist./1e9, vre
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
            uvdist./1e9, vim
    end

    @series begin
        xguide --> "uv-distance (Gλ)"
        subplot := 3
        yguide := "LR (Jy)"
        vre = real.(getindex.(coh, 2, 1))
        vim = imag.(getindex.(coh, 2, 1))
        err = getindex.(noise, 1, 1)
        @series begin
            seriestype := :scatter
            alpha := 0.5
            yerr := err
            linecolor := nothing
            label := "Real"
            uvdist./1e9, vre
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
            uvdist./1e9, vim
    end
    @series begin
        xguide --> "uv-distance (Gλ)"
        subplot := 4
        yguide := "LL (Jy)"
        vre = real.(getindex.(coh, 2, 2))
        vim = imag.(getindex.(coh, 2, 2))
        err = getindex.(noise, 1, 1)
        @series begin
            seriestype := :scatter
            alpha := 0.5
            yerr := err
            linecolor := nothing
            label := "Real"
            uvdist./1e9, vre
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
            uvdist./1e9, vim
    end

end

@recipe function f(dvis::EHTObservationTable{A};) where {A<:EHTVisibilityAmplitudeDatum}
    xguide --> "uv-distance (Gλ)"
    yguide --> "|V| (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :U)/1e9
    v = getdata(dvis, :V)/1e9
    uvdist = hypot.(u,v)
    amp = dvis[:measurement]
    noise = getdata(dvis, :noise)
    #add data noisebars
    seriestype --> :scatter
    alpha --> 0.5
    yerr := noise
    linecolor --> nothing
    label --> "Data"
    uvdist, amp
end

@recipe function f(acc::AbstractArrayConfiguration)
    xguide --> "u (Gλ)"
    yguide --> "v (Gλ)"
    markershape --> :circle

    u, v = getuv(acc)
    #add data noisebars
    seriestype --> :scatter
    linecolor --> nothing
    aspect_ratio --> :equal
    label -->"Data"
    title --> "Frequency: $(first(acc.data.F)/1e9) GHz"
    vcat(u/1e9,-u/1e9), vcat(v/1e9,-v/1e9)
end

@recipe function f(m::AbstractModel, dvis::EHTObservationTable{A}; datamarker=:circle, datacolor=:grey) where {A<:EHTVisibilityAmplitudeDatum}
    xguide --> "uv-distance (Gλ)"
    yguide --> "V (Jy)"
    markershape --> :diamond

    u = getdata(dvis, :U)
    v = getdata(dvis, :V)
    uvdist = hypot.(u,v)
    amp = dvis[:measurement]
    noise = getdata(dvis, :error)
    #add data noisebars
    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := datacolor
        alpha := 0.5
        yerr := noise
        linecolor := nothing
        label := "Data"
        uvdist./1e9, amp
    end

    seriestype-->:scatter
    amod = abs.(visibilitymap(m, arrayconfig(dvis)))
    labels --> "Model"
    uvdist./1e9, amod
end

export uvarea
function uvarea(d::EHTClosurePhaseDatum)
    u = map(x->x.U, d.baseline)
    v = map(x->x.V, d.baseline)
    a = hypot(u[1]-u[2], v[1]-v[2])
    b = hypot(u[2]-u[3], v[2]-v[3])
    c = hypot(u[3]-u[1], v[3]-v[1])
    heron(a,b,c)
end

function heron(a,b,c)
    s = 0.5*(a+b+c)
    return sqrt(s*(s-a)*(s-b)*(s-c))
end

function uvarea(d::EHTLogClosureAmplitudeDatum)
    u = map(x->x.U, d.baseline)
    v = map(x->x.V, d.baseline)
    a = hypot(u[1]-u[2], v[1]-v[2])
    b = hypot(u[2]-u[3], v[2]-v[3])
    c = hypot(u[3]-u[4], v[3]-v[4])
    d = hypot(u[4]-u[1], v[4]-v[1])
    h = hypot(u[1]-u[3], v[1]-v[3])
    return heron(a,b,h)+heron(c,d,h)
end

@recipe function f(dlca::EHTObservationTable{A}) where {A<:EHTLogClosureAmplitudeDatum}
    xguide --> "√(quadrangle area) (λ)"
    yguide --> "Log Clos. Amp."
    markershape --> :diamond
    area = sqrt.(uvarea.(datatable(dlca)))
    phase = measurement(dlca)
    err = noise(dlca)
    #add data noisebars
    seriestype --> :scatter
    alpha := 0.5
    yerr := err
    label --> "Data"
    linecolor --> nothing
    area, phase
end


@recipe function f(m::AbstractModel, dlca::EHTObservationTable{A}; datamarker=:circle, datacolor=:grey) where {A<:EHTLogClosureAmplitudeDatum}
    xguide --> "√(quadrangle area) (λ)"
    yguide --> "Log Clos. Amp."
    markershape --> :diamond
    area = sqrt.(uvarea.(datatable(dlca)))
    phase = measurement(dlca)
    noise = noise(dlca)
    #add data noisebars
    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := datacolor
        alpha := 0.5
        yerr := noise
        linecolor := nothing
        label := "Data"
        area, phase
    end

    seriestype-->:scatter
    amod = logclosure_amplitudes(m, dlca.config)
    labels --> "Model"
    area,amod
end

@recipe function f(dcp::EHTObservationTable{A}) where {A<:EHTClosurePhaseDatum}
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
    noise = getdata(dcp, :noise)
    seriestype := :scatter
    alpha --> 0.5
    yerr := noise
    linecolor --> nothing
    label --> "Data"
    area, phase
end

@recipe function f(m::AbstractModel, dcp::EHTObservationTable{A}; datamarker=:circle, datacolor=:grey) where {A<:EHTClosurePhaseDatum}
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
    noise = getdata(dcp, :noise)
    #add data noisebars
    @series begin
        seriestype := :scatter
        markershape := datamarker
        markercolor := datacolor
        alpha := 0.5
        yerr := noise
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
ndata(d::EHTObservationTable) = length(d)
ndata(d::EHTObservationTable{D}) where {D<:EHTVisibilityDatum} = 2*length(d)
ndata(d::EHTObservationTable{D}) where {D<:EHTCoherencyDatum} = 8*length(d)

@recipe function f(h::Residual)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: AbstractModel) ||
        !(typeof(h.args[2]) <: EHTObservationTable)
        noise("Residual should be given a model and data product.  Got: $(typeof(h.args))")
    end
    m, damp = h.args
    uvdist, res = residuals(m, damp)
    c2 = chi2(m, damp)
    # title-->"Norm. Residuals"
    legend-->nothing

    if damp isa EHTObservationTable{<:EHTCoherencyDatum}
        layout := (2,2)
        res2 = reinterpret(reshape, Float64, res)'
        @series begin
            yguide := "RR"
            subplot := 1
            seriestype := :scatter
            alpha := 0.5
            linecolor := nothing
            title --> @sprintf "χ² = %.2f" sum(abs2, filter(!isnan, @view res2[:,1:2]))/(2*size(res2,1))
            uvdist./1e9, res2[:,1:2]
        end
        @series begin
            xguide --> "uv-distance (Gλ)"
            yguide := "LR"
            subplot := 3
            seriestype := :scatter
            alpha := 0.5
            linecolor := nothing
            title --> @sprintf "χ² = %.2f" sum(abs2, filter(!isnan, @view res2[:,3:4]))/(2*size(res2,1))
            uvdist./1e9, res2[:,3:4]
        end
        @series begin
            yguide := "RL"
            subplot := 2
            seriestype := :scatter
            alpha := 0.5
            linecolor := nothing
            title --> @sprintf "χ² = %.2f" sum(abs2, filter(!isnan, @view res2[:,5:6]))/(2*size(res2,1))
            uvdist./1e9, res2[:,5:6]
        end
        @series begin
            yguide := "LL"
            subplot := 4
            seriestype := :scatter
            alpha := 0.5
            linecolor := nothing
            title --> @sprintf "χ² = %.2f" sum(abs2, filter(!isnan, @view res2[:,7:8]))/(2*size(res2,1))
            uvdist./1e9, res2[:,7:8]
        end
    else
        xguide --> "uv-distance (Gλ)"
        yguide --> "Normalized Residual"
        markershape --> :circle
        linecolor --> nothing
        legend --> false

        title --> @sprintf "<χ²> = %.2f" c2/ndata(damp)
        return uvdist./1e9, res
    end
end


function chi2(m, data::EHTObservationTable)
    return sum(x->abs2.(x), last(residuals(m, data)))
end

function chi2(m, data::EHTObservationTable{A}) where {A<:EHTCoherencyDatum}
    res = last(residuals(m, data))
    return mapreduce(+, 1:4) do i
        sum(abs2, filter(!isnan, getproperty(res, i)))
    end
end


function chi2(m, data::EHTObservationTable...)
    return mapreduce(d->chi2(m, d), +, data)
end


function residuals_data(vis, data::EHTObservationTable{A}) where {A<:EHTClosurePhaseDatum}
    phase = measurement(data)
    err = noise(data)

    mphase = closure_phases(vis, arrayconfig(data))
    res = @. abs(cis(phase) - cis(mphase))
    return EHTObservationTable{A}(res, err, arrayconfig(data))
end


function residuals_data(vis, dlca::EHTObservationTable{A}) where {A<:EHTLogClosureAmplitudeDatum}
    phase = measurement(dlca)
    err = noise(dlca)
    mphase = logclosure_amplitudes(vis, arrayconfig(dlca))
    res = (phase .- mphase)
    return EHTObservationTable{A}(res, err, arrayconfig(dlca))
end

function residuals_data(vis, damp::EHTObservationTable{A}) where {A<:EHTVisibilityAmplitudeDatum}
    mamp = abs.(vis)
    amp = measurement(damp)
    res = (amp - mamp)
    return EHTObservationTable{A}(res, noise(damp), arrayconfig(damp))
end

function residuals_data(vis, dvis::EHTObservationTable{A}) where {A<:EHTCoherencyDatum}
    coh = measurement(dvis)
    res = coh .- vis
    return EHTObservationTable{A}(res, noise(dvis), arrayconfig(dvis))
end

function residuals_data(mvis, dvis::EHTObservationTable{A}) where {A<:EHTVisibilityDatum}
    vis = measurement(dvis)
    res = (vis - mvis)./getdata(dvis, :noise)
    return hypot.(u, v), res
end
