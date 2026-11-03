export residualplot

"""
    residualplot(post::AbstractVLBIPosterior, p; figure = (;), kwargs...)

A figure of the normalized residuals of every data product of the posterior `post` at the
parameters `p`, one row per product, against the baseline length. `figure` is passed to
`Figure` and the other keywords to [`fieldplot`](@ref). Needs a Makie backend.
"""
function residualplot end
