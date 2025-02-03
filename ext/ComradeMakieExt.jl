module ComradeMakieExt

using Makie
using Comrade
using Printf
import Measurements

import Comrade: plotfields, axisfields, plotcaltable, plotaxis
import Comrade: baselineplot, baselineplot!

@doc"""
    baselineplot(data, bl, fieldx, fieldy; kwargs...)
    baselineplot(data, fieldx, fieldy; kwargs...)

Plots the baseline `bl` with `fieldx` on the x axis and 
`fieldy` on the y axis.


If `bl` is a `Colon` or bl is not specified all baselines are plotted.

# Arguments
 - `data` : The data to plot

 - `field1` and `field2` : The fields to plot. field1 - x axis, field2 - y axis 
    If `field1` or `field2` is a function it will apply to `datatable(obsdata)` to get the value.
    If field1 or field2 is a symbol, it will look for a predefined function:
     - :U - baseline u coordinate
     - :V - baseline v coordinate
     - :Ti - time
     - :Fr - frequency
     - :measure - measurement
     - :noise - noise
     - :amp - visibility amplitude
     - :phase - visibility phase
     - :uvdist - projected baseline length
     - :snr - signal to noise ratio
     - :res - normalized residual visibilities (only if obsdata contains the residuals)

# Specific Attributes
  - `error` : If true, error bars are plotted. **This assumes that the data is a `Measurement` type.**
  - `color` : The color of the markers
  - `colorim` : The color of the imaginary part of the data if it exists
  - `marker` : The marker type
  - `markersize` : The size of the markers
  - `alpha` : The transparency of the markers
"""
Makie.@recipe(BaselinePlot, data, bl, fieldx, fieldy) do scene
    Makie.Attributes(;
        color = Makie.wong_colors()[1],
        colorim= Makie.wong_colors()[2],
        marker = :x,
        markersize = 10.0,
        alpha = 1.0,
        error = false,
    )
end

function Makie.convert_arguments(::Type{<:BaselinePlot}, data::Comrade.EHTObservationTable, fieldx, fieldy)
    return data, Colon(), fieldx, fieldy
end

# This is convienence function to conver the field to a function
# if the field is a symbol, we see if it is an option and return that
convert_field(f) = f
convert_field(::Type{<:U}) = x->x.baseline.U
convert_field(::Type{<:V}) = x->x.baseline.V
convert_field(::Type{<:Ti}) = x->x.baseline.Ti
convert_field(::Type{<:Fr}) = x->x.baseline.Fr

function convert_field(field::Symbol)
    field == :U && return x->x.baseline.U
    field == :V && return x->x.baseline.V
    field == :Ti && return x->x.baseline.Ti
    field == :Fr && return x->x.baseline.Fr
    field == :snr && return x->abs.(Comrade.measurement(x)) ./ noise(x)
    field == :uvdist && return uvdist
    field == :amp && return x->abs.(measurement(x))
    field == :phase && return x->angle.(measurement(x))
    field == :res && return x->Comrade.measurement(x) ./ noise(x)
    field == :measurement && return x->Comrade.measurement(x)
    field == :noise && return x->noise(x)
    field == :measwnoise && return x->measwnoise(x)

    throw(ArgumentError("$field not supported please use one of"* 
                        "(:U, :V, :Ti, :Fr, :snr, :uvdist, :amp, :phase,\n"*
                        " :res, :measurement, :noise, :measwnoise)"))
end

function Makie.plot!(plot::BaselinePlot{<:Tuple{<:Comrade.EHTObservationTable, 
                                                <:Union{Tuple{Symbol, Symbol}, Colon}, 
                                                <:Any, <:Any}})
    @extract plot (data, bl, fieldx, fieldy)
    obl = @lift begin
        if $bl isa Colon
            return datatable($data)
        else
            datatable(select_baseline($data, $bl))
        end
    end
    @lift begin
        if isempty($obl)
            throw(ArgumentError("Baseline $($bl) does not exist"))
        end
    end

    fx = @lift(convert_field($fieldx))
    fy = @lift(convert_field($fieldy))

    x  = @lift($(fx).($obl))
    y  = @lift($(fy).($obl))

    lift(x, y, plot.error) do x, y, berr
        if eltype(y) <: Complex
            Makie.scatter!(plot, x, real.(y); 
                            color=plot.color, 
                            marker=plot.marker, 
                            markersize=plot.markersize, 
                            alpha=plot.alpha)
            Makie.scatter!(plot, x, imag.(y); 
                            color=plot.colorim, 
                            marker=plot.marker, 
                            markersize=plot.markersize, 
                            alpha=plot.alpha)

            if berr
                eltype(y) <: Complex{<:Measurements.Measurement} || throw(ArgumentError("Field does not have error, did use measwnoise?"))
                Makie.errorbars!(plot, x, real.(y), 
                                 color = plot.color,
                                 alpha = plot.alpha
                                )
                Makie.errorbars!(plot, x, imag.(y), 
                                 color = plot.colorim,
                                 alpha = plot.alpha
                                )

            end


        else
            Makie.scatter!(plot, x, y; 
                            color=plot.color, 
                            marker=plot.marker, 
                            markersize=plot.markersize, 
                            alpha=plot.alpha)
            if berr
                eltype(y) <: Measurements.Measurement || throw(ArgumentError("Field does not have error, did use measwnoise?"))
                Makie.errorbars!(plot, x, y, 
                                 color = plot.color,
                                 alpha = plot.alpha
                                )
            end
                            
        end
    end
    return plot
end



function frequencylabel(ν::Number)
    if ν > 1e6
        if ν > 1e12
            label = @sprintf("%.2f", ν / 1e12) * " THz"
        elseif ν > 1e9
            label = @sprintf("%.2f", ν / 1e9) * " GHz"
        else
            label = @sprintf("%.2f", ν / 1e6) * " MHz"
        end
    else
        label = @sprintf("%.2f", ν) * " Hz"
    end
    return label
end

_frequency(d::EHTObservationTable) = domain(d).Fr
_frequency(d::EHTObservationTable{<:Comrade.ClosureProducts}) = datatable(d).baseline.:(1).Fr
measname(d::EHTObservationTable{<:Comrade.EHTVisibilityDatum}) = "Visibility (Jy)"
measname(d::EHTObservationTable{<:Comrade.EHTClosurePhaseDatum}) = "Closure Phase (rad)"
measname(d::EHTObservationTable{<:Comrade.EHTLogClosureAmplitudeDatum}) = "Log Closure Amplitude"
measname(d::EHTObservationTable{<:Comrade.EHTCoherencyDatum}) = "Coherency (Jy)"
measname(d::EHTObservationTable{<:Comrade.EHTVisibilityAmplitudeDatum}) = "Visibility Amplitude (Jy)"

_defaultlabel(f::Symbol) = _defaultlabel(Val(f))
_defaultlabel(::Val{:U}) = "u (λ)"
_defaultlabel(::Val{:V}) = "v (λ)"
_defaultlabel(::Val{:Ti}) = "Time (UTC)"
_defaultlabel(::Val{:Fr}) = "Frequency (Hz)"
_defaultlabel(::Val{:amp}) = "Visibility Amplitude (Jy)"
_defaultlabel(::Val{:phase}) = "Visibility Phase (rad)"
_defaultlabel(::Val{:uvdist}) = "Projected Baseline Distance λ"
_defaultlabel(::Val{:snr}) = "SNR"
_defaultlabel(::Val{:res}) = "Normalized Residual Visibility"
_defaultlabel(::typeof(Comrade.uvdist)) = _defaultlabel(:uvdist)
_defaultlabel(f) = string(f)

@doc"""

    plotfields(obsdata::EHTObservationTable, field1, field2; 
                site1=nothing, site2=nothing, axis_kwargs=(;), legend_kwargs=(;), scatter_kwargs=(;))

Plots two data fields against each other.

# Arguments
 - `obsdata` : EHTObservationTable containing the data to plot (closure quantities not supported yet)

 - `field1` and `field2` : The fields to plot. field1 - x axis, field2 - y axis 
    If field1 or field2 is a function it will apply to `datatable(obsdata)` to get the value
    If field1 or field2 is a symbol, it will look for a predefined function:
     - :U - baseline u coordinate
     - :V - baseline v coordinate
     - :Ti - time
     - :Fr - frequency
     - :measure - measurement
     - :noise - noise
     - :amp - visibility amplitude
     - :phase - visibility phase
     - :uvdist - projected baseline length
     - :snr - signal to noise ratio
     - :res - normalized residual visibilities (only if obsdata contains the residuals)

 - `site1` and `site2` : Keywords for the sites forming the baseline being plotted, e.g. :ALMA, :APEX.
 - `axis_kwargs` : Keyword arguments for each subplot's Axis.
 - `legend_kwargs` : Keyword arguments passed to the figure legend.
 - `scatter_kwargs` : Keyword arguments passed to scatter! in each subplot.
"""
function plotfields(
    obsdata::Comrade.EHTObservationTable,
    field1,
    field2;
    legend = true,
    conjugate = true,
    axis_kwargs = (;),
    legend_kwargs = (;),
    scatter_kwargs = (;),
)

    if field1 == :measurement
        xlabel = measname(obsdata)
    else 
        xlabel = _defaultlabel(field1)
    end

    if field2 == :measurement
        ylabel = measname(obsdata)
    else 
        ylabel = _defaultlabel(field2)
    end
    
    axis_kwargs = (;
        axis_kwargs...,
        xlabel = xlabel,
        ylabel = ylabel,
    )

    fx = convert_field(field1)
    fy = convert_field(field2)
    dt = datatable(obsdata)
    x = fx.(dt)
    y = fy.(dt)
    Fr = _frequency(obsdata)
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
    obsdata::Comrade.EHTObservationTable,
    field1,
    field2,
    site1::Symbol,
    site2::Symbol;
    axis_kwargs = (;),
    legend_kwargs = (;),
    scatter_kwargs = (;),
)
    title = string(site1) * " - " * string(site2)
    blobs = select_baseline(obsdata, (site1, site2))
    return plotfields(
        blobs,
        field1,
        field2;
        axis_kwargs = (; axis_kwargs..., title = title),
        legend_kwargs = legend_kwargs,
        scatter_kwargs = scatter_kwargs,
    )
end

@doc"""

    axisfields(fig, obsdata::EHTObservationField, field1, field2;
                legend=true, conjugate=true, axis_kwargs=(;), legend_kwargs=(;), scatter_kwargs=(;))

Plots two data fields from `obsdata` against each other on `fig`, returns a Makie Axis which can 
be used to configure subplots.

# Arguments
- `fig`: The GridPosition i.e. `fig[i,j]` to plot the data on.
- `obsdata` : EHTObservationTable containing the data to plot (closure quantities not supported yet)

- `field1` and `field2` : The fields to plot. field1 - x axis, field2 - y axis 
    If field1 or field2 is a function it will apply to `datatable(obsdata)` to get the value
    If field1 or field2 is a symbol, it will look for a predefined function:
     - :U - baseline u coordinate
     - :V - baseline v coordinate
     - :Ti - time
     - :Fr - frequency
     - :measure - measurement
     - :noise - noise
     - :amp - visibility amplitude
     - :phase - visibility phase
     - :uvdist - projected baseline length
     - :snr - signal to noise ratio
     - :res - normalized residual visibilities (only if obsdata contains the residuals)

# Keyword Arguments
 - `legend` : If true, legend is shown. If false, legend is hidden.
 - `conjugate` : Only relevant if plotting (u,v) coverage. If true, data is conjugated. If false, data is plotted as is. 
 - `site1` and `site2` : Keywords for the sites forming the baseline being plotted, e.g. :ALMA, :APEX.
 - `axis_kwargs` : Keyword arguments for each subplot's Axis.
 - `legend_kwargs` : Keyword arguments passed to the figure legend.
 - `scatter_kwargs` : Keyword arguments passed to scatter! in each subplot.
"""
function axisfields(
    fig::GridPosition,
    obsdata::Comrade.EHTObservationTable,
    field1,
    field2;
    legend = true,
    conjugate = true,
    axis_kwargs = (;),
    legend_kwargs = (;),
    scatter_kwargs = (;),
)
    if field1 == :measurement
        xlabel = measname(obsdata)
    else 
        xlabel = _defaultlabel(field1)
    end

    if field2 == :measurement
        ylabel = measname(obsdata)
    else 
        ylabel = _defaultlabel(field2)
    end

    axis_kwargs = (;
        axis_kwargs...,
        xlabel = xlabel,
        ylabel = ylabel,
    )

    fx = convert_field(field1)
    fy = convert_field(field2)
    dt = datatable(obsdata)
    x = fx.(dt)
    y = fy.(dt)
    Fr = _frequency(obsdata)

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
    obsdata::Comrade.EHTObservationTable,
    field1,
    field2,
    site1::Symbol,
    site2::Symbol;
    legend = true,
    axis_kwargs = (;),
    legend_kwargs = (;),
    scatter_kwargs = (;),
)
    title = string(site1) * " - " * string(site2)
    blobs = select_baseline(obsdata, (site1, site2))
    return axisfields(
        fig,
        blobs,
        field1,
        field2;
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
        νind = findall(==(ν), Fr)
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

@doc"""
    plotcaltable(gt...; width=150, height=125, layout=nothing, markers=nothing, labels=nothing, 
                        axis_kwargs=(;), legend_kwargs=(;), figure_kwargs=(;), scatter_kwargs=(;))

Automatically generate a grid of subplots plotting the Comrade.CalTable information.
Each subplot corresponds to a different station in the array.

## Argments

 - `gt` : The CalTables to plot.

## Keyword Arguments
 - `width` : Subplot width
 - `height` : Subplot height
 - `layout` : Subplot layout (optional). 
 - `axis_kwargs` : Keyword arguments for each subplot's Axis.
 - `legend_kwargs` : Keyword arguments passed to the figure legend.
 - `figure_kwargs` : Keyword arguments passed to Figure().
 - `scatter_kwargs` : Keyword arguments passed to scatter! in each subplot.
"""
function plotcaltable(
    gt::Comrade.CalTable...;
    width = 150,
    height = 125,
    layout = nothing,
    markers = nothing,
    labels = nothing,
    axis_kwargs = (;),
    legend_kwargs = (;),
    figure_kwargs = (;),
    scatter_kwargs = (;),
)
    if gt isa Comrade.CalTable
        gt = (gt,)
    end
    sitelist = sites(argmax(x->length(sites(x)), gt))
    tup = maximum(maximum.(x->x[:Ti].t0, gt))
    tlo = minimum(minimum.(x->x[:Ti].t0, gt))

    lims = get(axis_kwargs, :limits, ((tlo, tup), nothing))
    lims = isnothing(lims[1]) ? ((tlo, tup), lims[2]) : lims
    ylabel = get(axis_kwargs, :ylabel, nothing)

    axis_kwargs = Base.structdiff(axis_kwargs, NamedTuple{(:limits, :ylabel)})

    if isnothing(layout)
        collen = 4
        rowlen = ceil(Int, length(sitelist) / collen)
    else
        (rowlen, collen) = layout
    end

    size = (width * (collen), height * (rowlen) + 50) # height + 50 is for the UTC label
    fig = Figure(; size, figure_kwargs...)
    axs = map(Iterators.Filter(x->collen*(x[1]-1) + x[2] <= length(sitelist), Iterators.product(1:rowlen, 1:collen))) do (i,j)
        n = collen*(i-1) + j
        sitelist[n] => Axis(fig[i, j];
             title=string(sitelist[n]),
             width=width,
             height=height,
             limits = lims,
             axis_kwargs...)
    end

    markers = isnothing(markers) ? collect(keys(Makie.default_marker_map())) : markers
    labels = isnothing(labels) ? repeat(" ", length(sitelist)) : labels

    for (i, (site, ax)) in pairs(axs)
        for (k, gi) in pairs(gt)
            νlist = unique(gi.Fr)
            for (j, ν) in pairs(νlist)
                νind = findall(==(ν), gi.Fr)
                x = getproperty.(gi.Ti, :t0)[νind]
                y = getproperty(gi, site)[νind]

                if eltype(y) >: Float64
                    scatter!(
                        ax,
                        x,
                        y,
                        label = string(labels[k], " ", frequencylabel(round(ν.central, digits = 2))),
                        marker = markers[j],
                        scatter_kwargs...,
                    )
                else
                    missingind = findall(x -> typeof(x) == Missing, y)
                    y[missingind] .= NaN
                    yval = getproperty.(y, :val)
                    yerr = getproperty.(y, :err)

                    errorbars!(ax, x, yval, yerr, scatter_kwargs...)
                    scatter!(
                        ax,
                        x,
                        yval,
                        label = string(labels[i], " ", frequencylabel(round(ν.central, digits = 2))),
                        marker = markers[j],
                        scatter_kwargs...,
                    )
                end
            end
        end


    end

    Legend(fig[1, collen+1], axs[end][2], legend_kwargs...)
    axs[rowlen][2].xlabel = "Time (UTC)"
    if !isnothing(ylabel)
        axs[rowlen][2].ylabel = ylabel
    end
    resize_to_layout!(fig)
    return fig
end

end