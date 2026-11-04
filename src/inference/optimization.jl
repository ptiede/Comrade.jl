export estimate_mode

"""
    estimate_mode(post::AbstractVLBIPosterior, opt, args...; initial_params=nothing, kwargs...)

The mode of the posterior `post`, found with the optimizer `opt`. The `args/kwargs` are
forwarded to the optimization package.

!!! warning
    This function won't have any methods until the optimization package is loaded.
    We recommend loading the `Optimization.jl` package.
"""
function estimate_mode end
