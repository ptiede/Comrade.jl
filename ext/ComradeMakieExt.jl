module ComradeMakieExt

using Makie
using Comrade
using Printf
using Accessors: @reset, @set
import Measurements

import Comrade: plotfields, plotfields!, axisfields, plotcaltable, plotaxis
import Comrade: baselineplot, baselineplot!

const MakieGrid = Union{Makie.GridPosition, Makie.GridSubposition}

@doc"""
    baselineplot(data, bl, fieldx, fieldy; kwargs...)
    baselineplot(data, fieldx, fieldy; kwargs...)

Plots the baseline `bl` with `fieldx` on the x axis and 
`fieldy` on the y axis.


If `bl` is a `Colon` or bl is not specified all baselines are plotted.

# Arguments
 - `data` : The data to plot

 - `field1` and `field2` : The fields to plot. field1 - x axis, field2 - y axis 
    If `field1` or `field2` is a function, baselineplot will broadcast them over `datatable(obsdata)` 
    to get the x and y values as the outputs from `field1.(datatable(obsdata))` and `field2.(datatable(obsdata))`.
    If field1 or field2 is a symbol, it will look for a predefined function:
     - :U - baseline u coordinate
     - :V - baseline v coordinate
     - :Ti - time
     - :Fr - frequency
     - :measurement - measurement
     - :measwnoise - measurement and noise
     - :noise - noise
     - :amp - visibility amplitude
     - :phase - visibility phase
     - :uvdist - projected baseline length
     - :snr - signal to noise ratio
     - :res - normalized residual visibilities (only if obsdata contains the residuals)


# Examples
```julia-repl
julia> baselineplot(log_closure_amplitudes, Comrade.uvdist, Comrade.measwnoise, error=true)
```
"""
Makie.@recipe BaselinePlot (data::Comrade.AbstractObservationTable, bl::Any, fieldx::Any, fieldy::Any) begin
    "The color of the marker"
    color = @inherit markercolor
    "The color to use for the imaginary component if fieldy returns a complex number. Note this does not cycle with the color attribute"
    colorim = @inherit linecolor
    "The marker type to use"
    marker = @inherit marker
    "Sets the size of the markers"
    markersize = @inherit markersize
    "Sets the color around the markers"
    strokecolor = @inherit markerstrokecolor
    "Sets the width of the color around the markers"
    strokewidth = @inherit markerstrokewidth
    "Deprecated will automatically plot errors if the field is a `Measurements.Measurement`"
    error = true
    "Width of the error bar"
    linewidth = @inherit linewidth
    "Cap of the error bar"
    linecap = @inherit linecap
    "Width of the error bar cap"
    linecapwidth = 0
    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_colormap_attributes()...
    cycle = [:color]
end

function Makie.convert_arguments(::Type{<:BaselinePlot}, data::Comrade.EHTObservationTable, fieldx, fieldy)
    return data, Colon(), fieldx, fieldy
end

# This is convienence function to conver the field to a function
# if the field is a symbol, we see if it is an option and return that
convert_field(f) = f
convert_field(::Type{<:U}) = _U
convert_field(::Type{<:V}) = _V
convert_field(::Type{<:Ti}) = x -> x.baseline.Ti
convert_field(::Type{<:Fr}) = x -> x.baseline.Fr

_U(x) = x.baseline.U
_V(x) = x.baseline.V

function convert_field(field::Symbol)
    field == :U && return _U
    field == :V && return _V
    field == :Ti && return x -> x.baseline.Ti
    field == :Fr && return x -> x.baseline.Fr
    field == :snr && return x -> abs.(Comrade.measurement(x)) ./ noise(x)
    field == :uvdist && return uvdist
    field == :amp && return x -> abs.(measurement(x))
    field == :phase && return x -> angle.(measurement(x))
    field == :res && return x -> Comrade.measurement(x) ./ noise(x)
    field == :measurement && return x -> Comrade.measurement(x)
    field == :noise && return x -> noise(x)
    field == :measwnoise && return x -> measwnoise(x)

    throw(
        ArgumentError(
            "$field not supported please use one of" *
                "(:U, :V, :Ti, :Fr, :snr, :uvdist, :amp, :phase,\n" *
                " :res, :measurement, :noise, :measwnoise)"
        )
    )
end

function ismeasure(::AbstractArray{T}) where {T<:Union{Complex{<:Measurements.Measurement}, Measurements.Measurement}}
    return true
end

ismeasure(x) = false

function Makie.plot!(
        plot::BaselinePlot{
            <:Tuple{
                <:Comrade.EHTObservationTable,
                <:Union{Tuple{Symbol, Symbol}, Colon},
                <:Any, <:Any,
            },
        }
    )

    map!(plot, [:bl, :data], :obl) do bl, data
        if bl isa Colon
            obs = datatable(data)
        else
            obs = datatable(select_baseline(data, bl))
        end
        isempty(obs) && throw(ArgumentError("No data for baseline $bl"))
        return obs
    end

    map!(plot, [:fieldx, :fieldy, :obl], [:x, :y]) do fieldx, fieldy, obl
        x, y = convert_field(fieldx).(obl), convert_field(fieldy).(obl)
        return x, y
    end



    lift(plot.x, plot.y) do x, y
        if eltype(y) <: Complex

            if ismeasure(y) || ismeasure(x)
                Makie.errorbars!(
                    plot, plot.attributes, x, real.(y),
                    color = plot.color,
                    linewidth = plot.linewidth,
                    linecap = plot.linecap,
                    whiskerwidth = plot.linecapwidth,

                )
                Makie.errorbars!(
                    plot, plot.attributes, x, imag.(y),
                    color = plot.colorim,
                    linewidth = plot.linewidth,
                    linecap = plot.linecap,
                    whiskerwidth = plot.linecapwidth,
                )
            end


            Makie.scatter!(
                plot, plot.attributes, x, real.(y);
                color = plot.color,
                marker = plot.marker,
                markersize = plot.markersize,
                strokecolor = plot.strokecolor,
                strokewidth = plot.strokewidth,
            )
            Makie.scatter!(
                plot, plot.attributes, x, imag.(y);
                color = plot.colorim,
                marker = plot.marker,
                markersize = plot.markersize,
                strokecolor = plot.strokecolor,
                strokewidth = plot.strokewidth,
            )



        else
            if ismeasure(y) || ismeasure(x)
                Makie.errorbars!(
                    plot, plot.attributes, x, y,
                    color = plot.color,
                    linewidth = plot.linewidth,
                    linecap = plot.linecap,
                    whiskerwidth = plot.linecapwidth,

                )
            end

            Makie.scatter!(
                plot, plot.attributes, x, y;
                color = plot.color,
                marker = plot.marker,
                markersize = plot.markersize,
                strokecolor = plot.strokecolor,
                strokewidth = plot.strokewidth,
            )

        end
    end
    return plot
end


function frequencylabel(ν::Number)
    if ν > 1.0e6
        if ν > 1.0e12
            label = @sprintf("%.2f", ν / 1.0e12) * " THz"
        elseif ν > 1.0e9
            label = @sprintf("%.2f", ν / 1.0e9) * " GHz"
        else
            label = @sprintf("%.2f", ν / 1.0e6) * " MHz"
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
_defaultlabel(::typeof(_U)) = "u (λ)"
_defaultlabel(::typeof(_V)) = "v (λ)"
_defaultlabel(::Type{<:U}) = "u (λ)"
_defaultlabel(::Type{<:V}) = "v (λ)"

_defaultlabel(::Val{:Ti}) = "Time (UTC)"
_defaultlabel(::Type{<:Ti}) = "Time (UTC)"
_defaultlabel(::Val{:Fr}) = "Frequency (Hz)"
_defaultlabel(::Type{<:Fr}) = "Frequency (Hz)"
_defaultlabel(::Val{:amp}) = "Visibility Amplitude (Jy)"
_defaultlabel(::Val{:phase}) = "Visibility Phase (rad)"
_defaultlabel(::Val{:uvdist}) = "Projected Baseline Distance (λ)"
_defaultlabel(::Val{:snr}) = "SNR"
_defaultlabel(::Val{:res}) = "Normalized Residual"
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
        legend_kwargs = (;merge=true, unique=true),
        scatter_kwargs = (;),
        figure_kwargs = (;)
    )


    fig = Figure(;figure_kwargs...)
    ax = plotfields!(
        fig[1, 1],
        obsdata, field1, field2;
        legend, conjugate,
        axis_kwargs, 
        scatter_kwargs,
        legend_kwargs
    )
    return Makie.FigureAxisPlot(fig, ax.axis, ax.plot)
end

function axis_labels(field1, field2, axis_kwargs)
    hasproperty(axis_kwargs, :xlabel) || (axis_kwargs = merge(axis_kwargs, (xlabel = _defaultlabel(field1),)))
    hasproperty(axis_kwargs, :ylabel) || (axis_kwargs = merge(axis_kwargs, (ylabel = _defaultlabel(field2),)))
    return axis_kwargs
end

function maybeset(nt, key, val)
    hasproperty(nt, key) && return nt 
    return merge(nt, NamedTuple{(key,)}((val,)))
end


function plotfields!(f::MakieGrid,
        obsdata::Comrade.EHTObservationTable,
        field1,
        field2;
        legend = true,
        conjugate = true,
        axis_kwargs = (;),
        legend_kwargs = (;merge=true, unique=true),
        scatter_kwargs = (;),
)


    ax = plotaxis(
        f,
        obsdata,
        field1, 
        field2;
        legend,
        conjugate,
        axis_kwargs,
        scatter_kwargs,
        legend_kwargs,
    )
    return ax
end



function plotfields(
        obsdata::Comrade.EHTObservationTable,
        bl::NTuple{2, Union{String, Symbol}},
        field1,
        field2;
        conjugate = false,
        legend = true,
        axis_kwargs = (;),
        legend_kwargs = (;merge=true, unique=true),
        scatter_kwargs = (;),
    )
    blobs = select_baseline(obsdata, bl)
    title = string(bl[1]) * " - " * string(bl[2])
    return plotfields(
        blobs,
        field1,
        field2;
        legend, 
        conjugate,
        axis_kwargs = maybeset(axis_kwargs, :title, title),
        legend_kwargs = legend_kwargs,
        scatter_kwargs = scatter_kwargs,
    )
end

function plotfields!(f::MakieGrid,
        obsdata::Comrade.EHTObservationTable,
        bl::NTuple{2, Union{String, Symbol}},
        field1,
        field2;
        legend = true,
        conjugate = true,
        axis_kwargs = (;),
        legend_kwargs = (;merge=true, unique=true),
        scatter_kwargs = (;),
)

    blobs = select_baseline(obsdata, bl)
    title = string(bl[1]) * " - " * string(bl[2])

    ax = plotaxis(
        f,
        blobs,
        field1, 
        field2;
        legend,
        conjugate,
        axis_kwargs=maybeset(axis_kwargs, :title, title),
        scatter_kwargs,
        legend_kwargs,
    )
    return ax
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
 - `axis_kwargs` : Keyword arguments for each subplot's Axis.
 - `legend_kwargs` : Keyword arguments passed to the figure legend.
 - `scatter_kwargs` : Keyword arguments passed to scatter! in each subplot.
"""
axisfields(args...; kwargs...) = plotfields!(args...; kwargs...)

getfr(x) = x.Fr
getfr(x::NTuple) = getfr(x[1])

_ismeas(field) = field in (:measurement, measurement, :measwnoise, measwnoise)

function plotaxis(
        fig::MakieGrid,
        obsdata, field1, field2;
        conjugate = false,
        legend = true,
        axis_kwargs = (;),
        scatter_kwargs = (;),
        legend_kwargs = (;merge=true, unique=true),
    )

    fx = (_ismeas(field1) ? measname(obsdata) : field1)
    fy = (_ismeas(field2) ? measname(obsdata) : field2)

    if field2 == :res
        title = "χ² = $(round(Comrade._chi2(obsdata)/Comrade.ndata(obsdata), digits=2))"
        axis_kwargs = maybeset(axis_kwargs, :title, title)
    end


    axis_kwargs = axis_labels(fx, fy, axis_kwargs)



    fx = convert_field(field1)
    fy = convert_field(field2)
    dt = datatable(obsdata)
    x = fx.(dt)
    y = fy.(dt)
    Fr = _frequency(obsdata)
    if conjugate # conjugating for (u,v) plotting
        if fx in (_U, _V) && fy in (_U, _V) # conjugating for (u,v) plotting
            x = [x; -x]
            y = [y; -y]
            Fr = [Fr; Fr]
            axis_kwargs = maybeset(axis_kwargs, :aspect, 1)
        end


    end

    if fx == _U
        axis_kwargs = maybeset(axis_kwargs, :xreversed, true)
    end


    νlist = unique(Fr) # multifrequency support

    ax = Axis(fig; axis_kwargs...)

    pl = nothing
    for ν in νlist
        obsν = filter(x->getfr(x.baseline)≈ν, obsdata)
        label = frequencylabel(ν)
        color = get(scatter_kwargs, :color, Makie.wong_colors()[1])
        pl = baselineplot!(ax, obsν, field1, field2; color, label, scatter_kwargs...)

        if conjugate && fx in (_U, _V) && fy in (_U, _V)
            baselineplot!(ax, obsν, x->-fx(x), x->-fy(x); color, label, scatter_kwargs...)
        end

    end
    if legend
        axislegend(ax; legend_kwargs...)
    end

    return Makie.AxisPlot(ax, pl)
end

"""
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
        tables...;
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

    gt = map(prepare_caltable, (tables...,))
    sitelist = sites(argmax(x -> length(sites(x)), gt))
    tup = maximum(maximum.(x -> x[:Ti].t0 + x[:Ti].dt / 2, gt))
    tlo = minimum(minimum.(x -> x[:Ti].t0 - x[:Ti].dt / 2, gt))

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
    axs = map(Iterators.Filter(x -> collen * (x[1] - 1) + x[2] <= length(sitelist), Iterators.product(1:rowlen, 1:collen))) do (i, j)
        n = collen * (i - 1) + j
        sitelist[n] => Axis(
            fig[i, j];
            title = string(sitelist[n]),
            width = width,
            height = height,
            limits = lims,
            axis_kwargs...
        )
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
                label = string(labels[k], " ", frequencylabel(round(ν.central, digits = 2)))
                marker = markers[j]
                if eltype(y) >: Float64
                    scatter!(
                        ax,
                        x,
                        y,
                        marker = marker,
                        label = label,
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
                        label = label,
                        marker = marker,
                        scatter_kwargs...,
                    )
                end
            end
        end


    end

    Legend(fig[1, collen + 1], axs[end][2], legend_kwargs...)
    axs[rowlen][2].xlabel = "Time (UTC)"
    if !isnothing(ylabel)
        axs[rowlen][2].ylabel = ylabel
    end
    resize_to_layout!(fig)
    return fig
end

prepare_caltable(gt::Comrade.CalTable) = gt
prepare_caltable(gt::Comrade.SiteArray) = caltable(gt)

end
