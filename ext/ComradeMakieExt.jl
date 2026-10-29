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
    If `field1` or `field2` is a function, baselineplot will broadcast them over `datarows(data)` 
    to get the x and y values as the outputs from `field1.(datarows(data))` and `field2.(datarows(data))`.
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
julia> baselineplot(log_closure_amplitudes(t), :uvdist, :measwnoise)
```
"""
Makie.@recipe BaselinePlot (data::AbstractVector, bl::Any, fieldx::Any, fieldy::Any) begin
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

# The tables are plotted through their rows, `datarows(t)`.
const _Table = Union{Comrade.VisTable, Comrade.ClosureTable}

Makie.convert_arguments(::Type{<:BaselinePlot}, rows::AbstractVector, fieldx, fieldy) = (rows, Colon(), fieldx, fieldy)
Makie.convert_arguments(::Type{<:BaselinePlot}, rows::AbstractVector, bl, fieldx, fieldy) = (rows, bl, fieldx, fieldy)
Makie.convert_arguments(P::Type{<:BaselinePlot}, t::_Table, args...) = Makie.convert_arguments(P, datarows(t), args...)

# The rows on the baseline `bl`, in either order.
function select_baseline(rows, bl)
    hasproperty(rows, :sites) || throw(ArgumentError("closures are not on one baseline"))
    s = Symbol.(bl)
    return filter(r -> r.sites == s || reverse(r.sites) == s, rows)
end

_measwnoise(m::Real, n) = Measurements.measurement(m, n)
_measwnoise(m::Complex, n) = complex(Measurements.measurement(real(m), n), Measurements.measurement(imag(m), n))
_measwnoise(m::AbstractMatrix, n::AbstractMatrix) = map(_measwnoise, m, n)

# This is convienence function to conver the field to a function
# if the field is a symbol, we see if it is an option and return that
convert_field(f) = f
convert_field(::Type{<:U}) = _U
convert_field(::Type{<:V}) = _V
convert_field(::Type{<:Ti}) = _Ti
convert_field(::Type{<:Fr}) = _Fr

_U(x) = x.U
_V(x) = x.V
_Ti(x) = x.Ti
_Fr(x) = x.Fr
_uvdist(x) = x.uvdist
_snr(x) = abs.(x.measurement) ./ x.noise
_amp(x) = abs.(x.measurement)
_phase(x) = angle.(x.measurement)
_res(x) = x.measurement ./ x.noise
_measurement(x) = x.measurement
_noise(x) = x.noise
_measwnoise(x) = _measwnoise(x.measurement, x.noise)

function convert_field(field::Symbol)
    field == :U && return _U
    field == :V && return _V
    field == :Ti && return _Ti
    field == :Fr && return _Fr
    field == :snr && return _snr
    field == :uvdist && return _uvdist
    field == :amp && return _amp
    field == :phase && return _phase
    field == :res && return _res
    field == :measurement && return _measurement
    field == :noise && return _noise
    field == :measwnoise && return _measwnoise

    throw(
        ArgumentError(
            "$field not supported please use one of" *
                "(:U, :V, :Ti, :Fr, :snr, :uvdist, :amp, :phase,\n" *
                " :res, :measurement, :noise, :measwnoise)"
        )
    )
end

function ismeasure(::AbstractArray{T}) where {T <: Union{Complex{<:Measurements.Measurement}, Measurements.Measurement}}
    return true
end

ismeasure(x) = false

function Makie.plot!(
        plot::BaselinePlot{
            <:Tuple{
                <:AbstractVector,
                <:Union{Tuple{Symbol, Symbol}, Colon},
                <:Any, <:Any,
            },
        }
    )

    map!(plot, [:bl, :data], :obl) do bl, data
        if bl isa Colon
            obs = data
        else
            obs = select_baseline(data, bl)
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
    ν >= 1.0e12 && return @sprintf("%.2f", ν / 1.0e12) * " THz"
    ν >= 1.0e9 && return @sprintf("%.2f", ν / 1.0e9) * " GHz"
    ν >= 1.0e6 && return @sprintf("%.2f", ν / 1.0e6) * " MHz"
    ν >= 1.0e3 && return @sprintf("%.2f", ν / 1.0e3) * " kHz"
    return @sprintf("%.2f", ν) * " Hz"
end

measname(t::Comrade.VisTable) = _measname(eltype(Comrade.measurement(t)))
_measname(::Type{<:Complex}) = "Visibility (Jy)"
_measname(::Type{<:Real}) = "Visibility Amplitude (Jy)"
_measname(::Type{<:AbstractMatrix}) = "Correlation (Jy)"
measname(ct::Comrade.ClosureTable) = ct.product isa Comrade.ClosurePhases ? "Closure Phase (rad)" : "Log Closure Amplitude"

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
_defaultlabel(f) = string(f)

@doc"""

    plotfields(obsdata, field1, field2; axis_kwargs=(;), legend_kwargs=(;), scatter_kwargs=(;))
    plotfields(obsdata, bl, field1, field2; kwargs...)

Plots two data fields against each other. A field whose values are matrices, such as the
measurements of a table of all correlations (`Correlations()`), gets one axis per correlation
on a 2×2 grid, and the figure is returned.

# Arguments
 - `obsdata` : the `VisTable` or `ClosureTable` to plot
 - `bl` : a baseline, e.g. `(:AA, :AX)`, to plot only its rows

 - `field1` and `field2` : The fields to plot. field1 - x axis, field2 - y axis 
    If field1 or field2 is a function it will apply to `datarows(obsdata)` to get the value
    If field1 or field2 is a symbol, it will look for a predefined function:
     - :U - baseline u coordinate
     - :V - baseline v coordinate
     - :Ti - time
     - :Fr - frequency
     - :measurement - measurement
     - :measwnoise - measurement with its noise as error bars
     - :noise - noise
     - :amp - visibility amplitude
     - :phase - visibility phase
     - :uvdist - projected baseline length
     - :snr - signal to noise ratio
     - :res - normalized residual visibilities (only if obsdata contains the residuals)

 - `axis_kwargs` : Keyword arguments for each subplot's Axis.
 - `legend_kwargs` : Keyword arguments passed to the figure legend.
 - `scatter_kwargs` : Keyword arguments passed to scatter! in each subplot.
"""

function plotfields(
        obsdata::_Table,
        field1,
        field2;
        legend = true,
        conjugate = true,
        axis_kwargs = (;),
        legend_kwargs = (; merge = true, unique = true),
        scatter_kwargs = (;),
        figure_kwargs = (;)
    )
    fig = Figure(; figure_kwargs...)
    ax = plotfields!(fig[1, 1], obsdata, field1, field2; legend, conjugate, axis_kwargs, scatter_kwargs, legend_kwargs)
    return _figureplot(fig, ax)
end

# A field whose values are matrices, such as the measurements of a table of all correlations, is
# plotted on a 2×2 grid of axes, one per correlation.
_matrixvalued(rows, field) = !isempty(rows) && convert_field(field)(first(rows)) isa AbstractMatrix
_figureplot(fig, ax::Makie.AxisPlot) = Makie.FigureAxisPlot(fig, ax.axis, ax.plot)
_figureplot(fig, ::AbstractMatrix) = fig

function axis_labels(field1, field2, axis_kwargs)
    hasproperty(axis_kwargs, :xlabel) || (axis_kwargs = merge(axis_kwargs, (xlabel = _defaultlabel(field1),)))
    hasproperty(axis_kwargs, :ylabel) || (axis_kwargs = merge(axis_kwargs, (ylabel = _defaultlabel(field2),)))
    return axis_kwargs
end

function maybeset(nt, key, val)
    hasproperty(nt, key) && return nt
    return merge(nt, NamedTuple{(key,)}((val,)))
end

"""
    plotfields!(f::Makie.GridPosition, obsdata, field1, field2;
                kwargs...)

Plots the two data fields in the grid position `f`. For more information about the arguments
see the [`plotfields`](@ref) function. Returns the `AxisPlot`, or for a field whose values are
matrices a 2×2 matrix of them, one per correlation.
"""
function plotfields!(f::MakieGrid, obsdata::_Table, field1, field2; kwargs...)
    return _plotrows!(f, datarows(obsdata), measname(obsdata), field1, field2; kwargs...)
end

function plotfields(
        obsdata::_Table,
        bl::NTuple{2, Union{String, Symbol}},
        field1,
        field2;
        conjugate = false,
        legend = true,
        axis_kwargs = (;),
        legend_kwargs = (; merge = true, unique = true),
        scatter_kwargs = (;),
        figure_kwargs = (;)
    )
    fig = Figure(; figure_kwargs...)
    ax = plotfields!(fig[1, 1], obsdata, bl, field1, field2; legend, conjugate, axis_kwargs, scatter_kwargs, legend_kwargs)
    return _figureplot(fig, ax)
end

function plotfields!(
        f::MakieGrid,
        obsdata::_Table,
        bl::NTuple{2, Union{String, Symbol}},
        field1,
        field2;
        axis_kwargs = (;),
        kwargs...
    )
    title = string(bl[1]) * " - " * string(bl[2])
    rows = select_baseline(datarows(obsdata), bl)
    return _plotrows!(f, rows, measname(obsdata), field1, field2; axis_kwargs = maybeset(axis_kwargs, :title, title), kwargs...)
end

function _plotrows!(f, rows, mname, field1, field2; kwargs...)
    _matrixvalued(rows, field2) || return plotaxis(f, rows, mname, field1, field2; kwargs...)
    g = GridLayout(f)
    return [plotaxis(g[i, j], rows, mname, field1, field2; entry = (i, j), kwargs...) for i in 1:2, j in 1:2]
end

@doc"""

    axisfields(fig, obsdata, field1, field2;
                legend=true, conjugate=true, axis_kwargs=(;), legend_kwargs=(;), scatter_kwargs=(;))

Plots two data fields from `obsdata` against each other on `fig`, returns a Makie Axis which can 
be used to configure subplots.

! warning
    This function is deprecated, use `plotfields!` instead.

"""
axisfields(args...; kwargs...) = plotfields!(args...; kwargs...)

# The entry `(i, j)` of a field whose values are matrices.
struct _Component{F}
    f::F
    i::Int
    j::Int
end
(c::_Component)(x) = c.f(x)[c.i, c.j]
_component(f, ::Nothing) = f
_component(f, (i, j)) = _Component(f, i, j)

# The reduced χ² of residual rows: Σ |r/σ|² over the measured components, per real datum.
function _chi2(rows, entry = nothing)
    r = Iterators.flatten(map(_component(_res, entry), rows))
    z = [x for x in r if !isnan(x)]
    return sum(abs2, z) / (eltype(z) <: Complex ? 2 * length(z) : length(z))
end

_ismeas(field) = field in (:measurement, :measwnoise)

# One axis of `field2` against `field1` over `rows`, whose measurements are `mname`; `entry`
# picks one correlation of matrix-valued fields.
function plotaxis(
        fig::MakieGrid,
        rows, mname, field1, field2;
        entry = nothing,
        conjugate = false,
        legend = true,
        axis_kwargs = (;),
        scatter_kwargs = (;),
        legend_kwargs = (; merge = true, unique = true),
    )
    lx = (_ismeas(field1) ? mname : field1)
    ly = (_ismeas(field2) ? mname : field2)

    if field2 == :res
        axis_kwargs = maybeset(axis_kwargs, :title, "χ² = $(round(_chi2(rows, entry), digits = 2))")
    end
    isnothing(entry) || (axis_kwargs = maybeset(axis_kwargs, :ylabel, string(_defaultlabel(ly), " ", entry)))
    axis_kwargs = axis_labels(lx, ly, axis_kwargs)

    fx = convert_field(field1)
    fy = _component(convert_field(field2), entry)
    uv = fx in (_U, _V) && fy in (_U, _V)
    conjugate && uv && (axis_kwargs = maybeset(axis_kwargs, :aspect, 1))
    fx == _U && (axis_kwargs = maybeset(axis_kwargs, :xreversed, true))

    ax = Axis(fig; axis_kwargs...)

    pl = nothing
    for ν in unique(rows.Fr)
        obsν = filter(x -> x.Fr ≈ ν, rows)
        label = frequencylabel(ν)
        color = get(scatter_kwargs, :color, Makie.wong_colors()[1])
        pl = baselineplot!(ax, obsν, fx, fy; color, label, scatter_kwargs...)
        if conjugate && uv
            baselineplot!(ax, obsν, x -> -fx(x), x -> -fy(x); color, label, scatter_kwargs...)
        end
    end
    legend && axislegend(ax; legend_kwargs...)

    return Makie.AxisPlot(ax, pl)
end

function Comrade.residual(post::Comrade.AbstractVLBIPosterior, p; figure_kwargs = (;), kwargs...)
    res = Comrade.residuals(post, p)
    fig = Figure(; size = (600, 350 * length(res)), figure_kwargs...)
    for (i, r) in enumerate(res)
        plotfields!(fig[i, 1], r, :uvdist, :res; kwargs...)
    end
    return fig
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
# `caltable` accepts both plain SiteArrays and hierarchical instrument samples
# (e.g. GaussMarkovSitePrior with fitted hyperparameters)
prepare_caltable(gt::Union{Comrade.SiteArray, Comrade.HierarchicalSiteSample}) = caltable(gt)

end
