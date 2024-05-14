"""
    comrade_opt(post::VLBIPosterior, opt, adtype=Optimization.NoAD(), args...; initial_params=nothing, kwargs...)

Optimize the posterior `post` using the `opt` optimizer. The `adtype` specifies the automatic differentiation.
The `args/kwargs` are forwarded to `the specific optimization package.

!!! warning
    This function won't have any methods until the optimization package is loaded.
    We recommend loading the `Optimization.jl` package.
"""
function comrade_opt end

function comrade_laplace end
