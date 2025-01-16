module ComradeCairoMakieExt

using CairoMakie


"""
Returns the indices of an EHTObservationTable associated with a given baseline.
"""
function getbaselineind(obsdata::EHTObservationTable, site1::Symbol, site2::Symbol)
    sitedata = datatable(obsdata).baseline.sites
    if site1 != site2
        baselineind = findall(x -> (site1 in x) && (site2 in x), sitedata)
        return baselineind
    else
        return println("site1 and site2 must be different to form a baseline.")
    end
end


"""
Gets the observation data associated with a field. General purpose function.

# Arguments

 - `field` : The field for which data is retrieved. 
    Current fields supported
    :U - baseline u coordinate
    :V - baseline v coordinate
    :Ti - time
    :Fr - frequency
    :amp - visibility amplitude
    :phase - visibility phase
    :uvdist - projected baseline length
    :snr - signal to noise ratio
    :res - normalized residual visibilities (only if obsdata contains the residuals)
"""
function getobsdatafield(obsdata::EHTObservationTable, field::Symbol)
    if field in (:amp, :phase, :uvdist, :snr, :res)
        if field == :amp # calculate visibility amplitudes
            vis = Comrade.measurement(obsdata)
            amps = abs.(vis)
            return amps
        elseif field == :phase # calculate visibility phases
            vis = Comrade.measurement(obsdata)
            phases = angle.(vis)
            return phases
        elseif field == :snr
            vis = Comrade.measurement(obsdata)
            sigma = Comrade.noise(obsdata)
            snr = abs.(vis) .* inv.(sigma)
            return snr
        elseif field == :uvdist # calculate uv Distance
            dt = datatable(obsdata)
            return uvdist.(dt)
        elseif field == :res
            vis = Comrade.measurement(obsdata)
            sigma = Comrade.noise(obsdata)
            res = vis .* inv.(sigma)
            return res
        end
    elseif field in (:U, :V, :Ti, :Fr)
        bls = datatable(obsdata).baseline
        return getproperty(bls, field)
    else
        println("field not supported")
    end
end


function frequencylabel(ν::Number)
    if ν > 1e6
        if ν > 1e12
            label = label = string(ν / 1e12) * " THz"
        elseif ν > 1e9
            label = label = string(ν / 1e9) * " GHz"
        else
            label = label = string(ν / 1e6) * " MHz"
        end
    else
        label = label = string(ν) * " Hz"
    end
    return label
end

function plotfields(
    obsdata::EHTObservationTable,
    field1::Symbol,
    field2::Symbol;
    legend = true,
    conjugate = true,
    axis_kwargs = (;),
    legend_kwargs = (;),
    scatter_kwargs = (;),
)

    labels = (;
        U = L"v $(\lambda)$",
        V = L"v $(\lambda)$",
        Ti = "Time (UTC)",
        Fr = "Frequency (Hz)",
        amp = "Visibility Amplitude (Jy)",
        phase = "Visibility Phase (rad)",
        uvdist = L"Projected Baseline Distance $(\lambda)$",
        snr = "SNR",
        res = "Normalized Residual Visibility",
    )

    axis_kwargs = (;
        axis_kwargs...,
        xlabel = getproperty(labels, field1),
        ylabel = getproperty(labels, field2),
    )

    x = getobsdatafield(obsdata, field1)
    y = getobsdatafield(obsdata, field2)
    Fr = domain(obsdata).Fr

    if conjugate == true # conjugating for (u,v) plotting
        if field1 in (:U, :V) && field2 in (:U, :V) # conjugating for (u,v) plotting 
            x = [x; -x]
            y = [y; -y]
            Fr = [Fr; Fr]
            axis_kwargs = (; axis_kwargs..., aspect = 1)
        end
    end

    νlist = unique(Fr) # multifrequency support

    fig = Figure()
    plotaxis(
        fig[1, 1],
        x,
        y,
        Fr,
        νlist;
        legend = legend,
        axis_kwargs = (; axis_kwargs...),
        scatter_kwargs = (; scatter_kwargs...),
        legend_kwargs = (; legend_kwargs...),
    )
    return fig
end

function plotfields(
    obsdata::EHTObservationTable,
    field1::Symbol,
    field2::Symbol,
    site1::Symbol,
    site2::Symbol;
    axis_kwargs = (;),
    legend_kwargs = (;),
    scatter_kwargs = (;),
)
    title = string(site1) * " - " * string(site2)
    siteind = getbaselineind(obsdata, site1, site2)
    return plotfields(
        obsdata[siteind],
        field1::Symbol,
        field2::Symbol;
        axis_kwargs = (; axis_kwargs..., title = title),
        legend_kwargs = legend_kwargs,
        scatter_kwargs = scatter_kwargs,
    )
end

function axisfields(
    fig::GridPosition,
    obsdata::EHTObservationTable,
    field1::Symbol,
    field2::Symbol;
    legend = true,
    conjugate = true,
    axis_kwargs = (;),
    legend_kwargs = (;),
    scatter_kwargs = (;),
)
    labels = (;
        U = L"v $(\lambda)$",
        V = L"v $(\lambda)$",
        Ti = "Time (UTC)",
        Fr = "Frequency (Hz)",
        amp = "Visibility Amplitude (Jy)",
        phase = "Visibility Phase (rad)",
        uvdist = L"Projected Baseline Distance $(\lambda)$",
        snr = "SNR",
        res = "Normalized Residual Visibility",
    )

    axis_kwargs = (;
        axis_kwargs...,
        xlabel = getproperty(labels, field1),
        ylabel = getproperty(labels, field2),
    )

    x = getobsdatafield(obsdata, field1)
    y = getobsdatafield(obsdata, field2)
    Fr = domain(obsdata).Fr

    if conjugate == true # conjugating for (u,v) plotting
        if field1 in (:U, :V) && field2 in (:U, :V) # conjugating for (u,v) plotting 
            x = [x; -x]
            y = [y; -y]
            Fr = [Fr; Fr]
            axis_kwargs = (; axis_kwargs..., aspect = 1)
        end
    end

    νlist = unique(Fr) # multifrequency support
    ax = plotaxis(
        fig,
        x,
        y,
        Fr,
        νlist;
        legend = legend,
        axis_kwargs = (; axis_kwargs...),
        scatter_kwargs = (; scatter_kwargs...),
        legend_kwargs = (; legend_kwargs...),
    )
    return ax
end

function axisfields(
    fig::GridPosition,
    obsdata::EHTObservationTable,
    field1::Symbol,
    field2::Symbol,
    site1::Symbol,
    site2::Symbol;
    legend = true,
    axis_kwargs = (;),
    legend_kwargs = (;),
    scatter_kwargs = (;),
)
    title = string(site1) * " - " * string(site2)
    siteind = getbaselineind(obsdata, site1, site2)
    return axisfields(
        fig,
        obsdata[siteind],
        field1::Symbol,
        field2::Symbol;
        legend = legend,
        axis_kwargs = (; axis_kwargs..., title = title),
        legend_kwargs = legend_kwargs,
        scatter_kwargs = scatter_kwargs,
    )
end

function plotaxis(
    fig::GridPosition,
    x::AbstractArray,
    y::AbstractArray,
    Fr::AbstractArray,
    νlist::AbstractArray;
    legend = true,
    axis_kwargs = (;),
    scatter_kwargs = (;),
    legend_kwargs = (;),
)
    ax = Axis(fig; axis_kwargs...)
    for ν in νlist
        νind = findall(x -> x == ν, Fr)
        label = frequencylabel(ν)

        if eltype(x) <: Complex
            scatter!(
                ax,
                real.(x[νind]),
                y[νind],
                label = label * " Real";
                scatter_kwargs...,
            )
            scatter!(
                ax,
                imag.(x[νind]),
                y[νind],
                label = label * " Imag";
                scatter_kwargs...,
            )
        elseif eltype(y) <: Complex
            scatter!(
                ax,
                x[νind],
                real.(y[νind]),
                label = label * " Real";
                scatter_kwargs...,
            )
            scatter!(
                ax,
                x[νind],
                imag.(y[νind]),
                label = label * " Imag";
                scatter_kwargs...,
            )
        else
            scatter!(ax, x[νind], y[νind], label = label; scatter_kwargs...)
        end
    end

    if legend == true
        axislegend(ax; legend_kwargs...)
    end
end

function plotcaltable(
    gt::Comrade.CalTable;
    width = 150,
    height = 125,
    layout = Nothing,
    axis_kwargs = (;),
    legend_kwargs = (;),
    figure_kwargs = (;),
    scatter_kwargs = (;),
)
    sitelist = sites(gt)

    if layout == Nothing
        collen = 4
        rowlen = ceil(length(sitelist) / collen)
    else
        (rowlen, collen) = layout
    end

    size = (width * (collen), height * (rowlen) + 50) # height + 50 is for the UTC label

    fig = Figure(; size, figure_kwargs...)
    for n in range(1, rowlen) # loop over every station
        for m in range(1, collen)
            ind = Int64(collen * (n - 1) + m)

            if ind <= length(sitelist)
                site = sites(gt)[ind]

                ax = Axis(
                    fig[Int64(n), Int64(m)],
                    title = string(site),
                    width = width,
                    height = height,
                    axis_kwargs...,
                )

                νlist = unique(gt.Fr)
                for ν in νlist
                    νind = findall(x -> x == ν, gt.Fr)
                    x = getproperty.(gt.Ti, :t0)[νind]
                    y = getproperty(gt, site)[νind]

                    if eltype(y) >: Float64
                        scatter!(
                            ax,
                            x,
                            y,
                            label = frequencylabel(round(ν.central, digits = 2)),
                            scatter_kwargs...,
                        )
                    else
                        missingind = findall(x -> typeof(x) == Missing, y)
                        y[missingind] .= NaN
                        yval = getproperty.(y, :val)
                        yerr = getproperty.(y, :err)

                        errorbars!(x, yval, yerr, scatter_kwargs...)
                        scatter!(
                            x,
                            yval,
                            label = frequencylabel(round(ν.central, digits = 2)),
                            scatter_kwargs...,
                        )
                    end
                end

                if n == 1 && m == 1
                    Legend(fig[1, Int64(collen + 1)], ax, legend_kwargs...)
                end

                if n == rowlen && m == 1
                    ax.xlabel = "Time (UTC)"
                end
            end
        end
    end
    resize_to_layout!(fig)
    fig
end

end
