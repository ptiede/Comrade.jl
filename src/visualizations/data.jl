export residuals, chi2

@recipe function f(dvis::EHTObservationTable{A};) where {A<:EHTVisibilityDatum}
    xguide --> "uv-distance (Gλ)"
    yguide --> "V (Jy)"
    markershape --> :circle

    u = dvis[:baseline].U
    v = dvis[:baseline].V
    uvdist = hypot.(u,v)
    vis = dvis[:measurement]
    noise = dvis[:noise]
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

    u = dvis[:baseline].U
    v = dvis[:baseline].V
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

    u = dvis[:baseline].U
    v = dvis[:baseline].V
    uvdist = hypot.(u,v)
    amp = dvis[:measurement]
    noise = getdata(dvis, :noise)
    #add data noisebars
    seriestype --> :scatter
    alpha --> 0.5
    yerr := noise
    linecolor --> nothing
    label --> "Data"
    uvdist./1e9, amp
end

@recipe function f(acc::AbstractArrayConfiguration)
    xguide --> "u (Gλ)"
    yguide --> "v (Gλ)"
    markershape --> :circle

    u = dvis[:baseline].U
    v = dvis[:baseline].V
    #add data noisebars
    seriestype --> :scatter
    linecolor --> nothing
    aspect_ratio --> :equal
    label -->"Data"
    title --> "Frequency: $(first(acc.data.F)/1e9) GHz"
    vcat(u/1e9,-u/1e9), vcat(v/1e9,-v/1e9)
end


export uvdist
uvdist(d) = hypot(d.U, d.V)

function uvdist(d::EHTClosurePhaseDatum)
    u = map(x->x.U, d.baseline)
    v = map(x->x.V, d.baseline)
    a = hypot(u[1]-u[2], v[1]-v[2])
    b = hypot(u[2]-u[3], v[2]-v[3])
    c = hypot(u[3]-u[1], v[3]-v[1])
    sqrt(heron(a,b,c))
end

function heron(a,b,c)
    s = 0.5*(a+b+c)
    return sqrt(s*(s-a)*(s-b)*(s-c))
end

function uvdist(d::EHTLogClosureAmplitudeDatum)
    u = map(x->x.U, d.baseline)
    v = map(x->x.V, d.baseline)
    a = hypot(u[1]-u[2], v[1]-v[2])
    b = hypot(u[2]-u[3], v[2]-v[3])
    c = hypot(u[3]-u[4], v[3]-v[4])
    d = hypot(u[4]-u[1], v[4]-v[1])
    h = hypot(u[1]-u[3], v[1]-v[3])
    return sqrt(heron(a,b,h)+heron(c,d,h))
end

@recipe function f(dlca::EHTObservationTable{A}) where {A<:EHTLogClosureAmplitudeDatum}
    xguide --> "√(convex quadrangle area) (λ)"
    yguide --> "Log Clos. Amp."
    markershape --> :diamond
    area = (uvdist.(datatable(dlca)))
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


@recipe function f(dcp::EHTObservationTable{A}) where {A<:EHTClosurePhaseDatum}
    xguide --> "√(triangle area) (λ)"
    yguide --> "Phase (rad)"
    markershape --> :circle
    area = (uvdist.(datatable(dcp)))
    phase = dcp[:measurement]
    noise = dcp[:noise]
    seriestype := :scatter
    alpha --> 0.5
    yerr := noise
    linecolor --> nothing
    label --> "Data"
    area, phase
end

@userplot Residual

export ndata
ndata(d::EHTObservationTable) = length(d)
ndata(d::EHTObservationTable{D}) where {D<:EHTVisibilityDatum} = 2*length(d)
ndata(d::EHTObservationTable{D}) where {D<:EHTCoherencyDatum} = 8*length(d)

@recipe function f(h::Residual)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: AbstractVLBIPosterior)
        noise("Residual should be given a posterior and parameters.  Got: $(typeof(h.args))")
    end
    post, p = h.args
    ress = residuals(post, p)
    c2s = chi2(post, p)
    # title-->"Norm. Residuals"
    legend-->nothing
    layout--> (length(ress), 1)
    size --> (600, 300*length(ress))

    for i in eachindex(ress)
        rest = ress[i]
        res = map(datatable(rest)) do d
            d.measurement./d.noise
        end
        c2 = c2s[i]
        uvdist = Comrade.uvdist.(datatable(rest))

        if res isa EHTObservationTable{<:EHTCoherencyDatum}
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
            @series begin
            xguide --> "uv-distance (Gλ)"
            T = datumtype(rest)
            ST = split("$T", "Datum{")[1] |> x->split(x, ".EHT")[end]
            yguide --> "Norm. Res. $ST"
            markershape --> :circle
            linecolor --> nothing
            legend --> false
            subplot := i
            title --> @sprintf "<χ²> = %.2f" c2/ndata(rest)
            if eltype(res) isa Complex
                res = reinterpret(reshape, Float64, res)'
            end
            return uvdist./1e9, res
            end
        end
    end
end


function chi2(post::AbstractVLBIPosterior, p)
    res = residuals(post, p)
    return map(_chi2, res)
end

function _chi2(res::EHTObservationTable)
    return sum(datatable(res)) do d
            return abs2(d.measurement/d.noise)
    end
end



function residuals_data(vis, data::EHTObservationTable{A}) where {A<:EHTClosurePhaseDatum}
    phase = measurement(data)
    err = noise(data)

    mphase = closure_phases(vis, arrayconfig(data))
    res = @. atan(sin(phase - mphase), cos(phase - mphase))
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
    res = (vis - mvis)
    return EHTObservationTable{A}(res, noise(dvis), arrayconfig(dvis))
end
