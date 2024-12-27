export getbaselineind, getobsdatafield, plotfields, axisfields
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

- `field` : Symbol. 

    The field to retrieve data. 
    
    Current fields supported:

    :U - baseline u coordinate
    
    :V - baseline v coordinate
    
    :Ti - time
    
    :Fr - frequency
    
    :amp - visibility amplitude
    
    :phase - visibility phase

    :uvdist - projected baseline length

    :snr - signal to noise ratio
"""
function getobsdatafield(obsdata::EHTObservationTable{}, field::Symbol)
    if field in (:amp, :phase, :uvdist, :snr)
        vis = obsdata.measurement
        if field == :amp # calculate visibility amplitudes
            amps = abs.(vis)
            return amps
        elseif field == :phase # calculate visibility phases
            phases = atan.(imag(vis), real(vis))
            return phases
        elseif field == :snr
            snr = abs.(vis) ./ obsdata.noise
            return snr
        elseif field == :uvdist # calculate uv Distance
            obsdomain = domain(obsdata)
            u = getfield.(getfield(obsdomain,:dims),:U)
            v = getfield.(getfield(obsdomain,:dims),:V)
            uvdist = hypot.(u,v)
            return uvdist
        end
    elseif field in (:U, :V, :Ti, :Fr)
        obsdomain = domain(obsdata)
        return getfield.(getfield(obsdomain,:dims),field)
    else
        println("field not supported")
    end
end


function frequencylabel(ν::Number)
    if ν > 1e6
        if ν > 1e12
            label = label=string(ν/1e12) * " THz"
        elseif ν > 1e9
            label = label=string(ν/1e9) * " GHz"
        else
            label = label=string(ν/1e6) * " MHz"
        end
    else
        label = label=string(ν) * " Hz"
    end
    return label
end


"""
Plots two data fields against each other, returns a Makie Figure.

# Arguments
- `obsdata` : EHTObservationTable containing the data to plot (closure quantities not supported yet)

- `field1` and `field2` : Symbol. 

    The fields to plot. field1 - x axis, field2 - y axis 
    
    Current fields supported:

    :U - baseline u coordinate
    
    :V - baseline v coordinate
    
    :Ti - time
    
    :Fr - frequency
    
    :amp - visibility amplitude
    
    :phase - visibility phase

    :uvdist - projected baseline length

    :snr - signal to noise ratio

- `site1` and `site2` : Symbol
    
    The sites forming the baseline to plot

- `axis_kwargs` : NamedTuple

    NamedTuple of the Axis keyword arguments

- `legend_kwargs` : NamedTuple

    NamedTuple of the axis legend keyword arguments

- `scatter_kwargs` : NamedTuple

    NamedTuple of the scatter plot keyword arguments
"""
function plotfields(obsdata::EHTObservationTable{}, field1::Symbol, field2::Symbol; axis_kwargs=(;), legend_kwargs=(;), scatter_kwargs=(;))

    labels = (; U=L"v $(\lambda)$", V=L"v $(\lambda)$", 
              Ti="Time (UTC)", Fr ="Frequency (Hz)", 
              amp="Visibility Amplitude (Jy)", phase="Visibility Phase (rad)", 
              uvdist=L"Projected Baseline Distance $(\lambda)$", snr="SNR")

    axis_kwargs = (; axis_kwargs..., xlabel=getfield(labels,field1), ylabel=getfield(labels,field2))

    x = getobsdatafield(obsdata, field1)
    y = getobsdatafield(obsdata, field2)
    Fr = domain(obsdata).Fr

    if field1 in (:U,:V) && field2 in (:U,:V) # conjugating for (u,v) plotting 
        x = [x;-x]
        y = [y;-y]
        Fr = [Fr;Fr]
        axis_kwargs = (;axis_kwargs..., aspect=1)
    end

    νlist = unique(Fr) # multifrequency support

    with_theme(theme_latexfonts()) do
        fig = Figure()#size=(500,400))
        ax1 = Axis(fig[1, 1]; axis_kwargs...)
        for ν in νlist
            νind = findall(x -> x == ν, Fr)
            label = frequencylabel(ν)
            scatter!(ax1, x[νind], y[νind], label=label; scatter_kwargs...)
        end

        axislegend(ax1; legend_kwargs...)

        return fig
    end
end

function plotfields(obsdata::EHTObservationTable{}, field1::Symbol, field2::Symbol, site1::Symbol, site2::Symbol; axis_kwargs=(;), legend_kwargs=(;), scatter_kwargs=(;))
    title = string(site1) * " - " * string(site2)
    siteind = getbaselineind(obsdata, site1, site2)
    return plotfields(obsdata[siteind], field1::Symbol, field2::Symbol; axis_kwargs=(;axis_kwargs..., title=title), legend_kwargs=legend_kwargs, scatter_kwargs=scatter_kwargs)
end

"""
Plots two data fields against each other, returns a Makie Axis which can be used to configure subplots.

# Arguments
- `obsdata` : EHTObservationTable containing the data to plot (closure quantities not supported yet)

- `field1` and `field2` : Symbol. 

    The fields to plot. field1 - x axis, field2 - y axis 
    
    Current fields supported:

    :U - baseline u coordinate
    
    :V - baseline v coordinate
    
    :Ti - time
    
    :Fr - frequency
    
    :amp - visibility amplitude
    
    :phase - visibility phase

    :uvdist - projected baseline length

    :snr - signal to noise ratio

- `site1` and `site2` : Symbol
    
    The sites forming the baseline to plot

- `axis_kwargs` : NamedTuple

    NamedTuple of the Axis keyword arguments

- `legend_kwargs` : NamedTuple

    NamedTuple of the axis legend keyword arguments

- `scatter_kwargs` : NamedTuple

    NamedTuple of the scatter plot keyword arguments
"""
function axisfields(fig::GridPosition, obsdata::EHTObservationTable{}, field1::Symbol, field2::Symbol; legend=true, axis_kwargs=(;), legend_kwargs=(;), scatter_kwargs=(;))
    set_theme!(theme_latexfonts())
    labels = (; U=L"v $(\lambda)$", V=L"v $(\lambda)$", 
              Ti="Time (UTC)", Fr ="Frequency (Hz)", 
              amp="Visibility Amplitude (Jy)", phase="Visibility Phase (rad)", 
              uvdist=L"Projected Baseline Distance $(\lambda)$", snr="SNR")

    axis_kwargs = (; axis_kwargs..., xlabel=getfield(labels,field1), ylabel=getfield(labels,field2))

    x = getobsdatafield(obsdata, field1)
    y = getobsdatafield(obsdata, field2)
    Fr = domain(obsdata).Fr

    if field1 in (:U,:V) && field2 in (:U,:V) # conjugating for (u,v) plotting 
        x = [x;-x]
        y = [y;-y]
        Fr = [Fr;Fr]
        axis_kwargs = (;axis_kwargs..., aspect=1)
    end

    νlist = unique(Fr) # multifrequency support

    ax = Axis(fig; axis_kwargs...)
    for ν in νlist
        νind = findall(x -> x == ν, Fr)
        label = frequencylabel(ν)
        scatter!(ax, x[νind], y[νind], label=label; scatter_kwargs...)
    end

    if legend == true
        axislegend(ax; legend_kwargs...)
    end

    return ax
end

