export plotfields

"""
Plots two data fields against each other.

# Arguments
 - `obsdata` : EHTObservationTable containing the data to plot (closure quantities not supported yet)

 - `field1` and `field2` : The fields to plot. field1 - x axis, field2 - y axis 
    Current fields supported:
    :U - baseline u coordinate
    :V - baseline v coordinate
    :Ti - time
    :Fr - frequency
    :amp - visibility amplitude
    :phase - visibility phase
    :uvdist - projected baseline length
    :snr - signal to noise ratio
    :res - normalized residual visibilities (only if obsdata contains the residuals)

 - `site1` and `site2` : Keywords for the sites forming the baseline being plotted, e.g. :ALMA, :APEX.
 - `axis_kwargs` : Keyword arguments for each subplot's Axis.
 - `legend_kwargs` : Keyword arguments passed to the figure legend.
 - `scatter_kwargs` : Keyword arguments passed to scatter! in each subplot.
"""
function plotfields end

"""
Plots two data fields against each other, returns a Makie Axis which can be used to configure subplots.

# Arguments
- `obsdata` : EHTObservationTable containing the data to plot (closure quantities not supported yet)

- `field1` and `field2` : The fields to plot. field1 - x axis, field2 - y axis 
    Current fields supported:
    :U - baseline u coordinate
    :V - baseline v coordinate
    :Ti - time
    :Fr - frequency
    :amp - visibility amplitude
    :phase - visibility phase
    :uvdist - projected baseline length
    :snr - signal to noise ratio
    :res - normalized residual visibilities (only if obsdata contains the residuals)

# Keyword Arguments
 - `legend` : If true, legend is shown. If false, legend is hidden.
 - `conjugate` : Only relevant if plotting (u,v) coverage. If true, data is conjugated. If false, data is plotted as is. 
 - `site1` and `site2` : Keywords for the sites forming the baseline being plotted, e.g. :ALMA, :APEX.
 - `axis_kwargs` : Keyword arguments for each subplot's Axis.
 - `legend_kwargs` : Keyword arguments passed to the figure legend.
 - `scatter_kwargs` : Keyword arguments passed to scatter! in each subplot.
"""
function axisfields end

"""
Automatically generate a grid of subplots plotting the Comrade.CalTable information.
Each subplot corresponds to a different station in the array.

## Argments

 - `gt` : The CalTable to plot.

## Keyword Arguments
 - `width` : Subplot width
 - `height` : Subplot height
 - `layout` : Subplot layout (optional). 
 - `axis_kwargs` : Keyword arguments for each subplot's Axis.
 - `legend_kwargs` : Keyword arguments passed to the figure legend.
 - `figure_kwargs` : Keyword arguments passed to Figure().
 - `scatter_kwargs` : Keyword arguments passed to scatter! in each subplot.
"""
function plotcaltable end