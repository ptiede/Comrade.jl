module ComradeMCMCDiagnosticToolsExt

using Comrade
using MCMCDiagnosticTools

"""
    ess(s::PosteriorSamples; kwargs...)

The effective sample size of every parameter of the draws × chains (or single chain) samples
`s`, nested like the output of [`Comrade.rmap`](@ref): a number for each scalar parameter and
an array of the parameter's shape for each array parameter. `kwargs` go to
`MCMCDiagnosticTools.ess`, e.g. `kind = :tail`.
"""
MCMCDiagnosticTools.ess(s::PosteriorSamples; kwargs...) = _diagnostic(ess, s; kwargs...)

"""
    rhat(s::PosteriorSamples; kwargs...)

The R-hat of every parameter of the draws × chains (or single chain) samples `s`, nested like
the output of [`Comrade.rmap`](@ref): a number for each scalar parameter and an array of the
parameter's shape for each array parameter. `kwargs` go to `MCMCDiagnosticTools.rhat`.
"""
MCMCDiagnosticTools.rhat(s::PosteriorSamples; kwargs...) = _diagnostic(rhat, s; kwargs...)

"""
    mcse(s::PosteriorSamples; kwargs...)

The Monte Carlo standard error of every parameter of the draws × chains (or single chain)
samples `s`, nested like the output of [`Comrade.rmap`](@ref): a number for each scalar
parameter and an array of the parameter's shape for each array parameter. `kwargs` go to
`MCMCDiagnosticTools.mcse`, e.g. `kind = std`.
"""
MCMCDiagnosticTools.mcse(s::PosteriorSamples; kwargs...) = _diagnostic(mcse, s; kwargs...)

function _diagnostic(f, s; kwargs...)
    ndims(s) <= 2 || throw(ArgumentError("diagnostics need draws × chains samples, not $(ndims(s)) dimensions"))
    return Comrade.rmap(x -> _per_parameter(f, x; kwargs...), s)
end

# MCMCDiagnosticTools takes a (draws, chains, parameters) array.
function _per_parameter(f, x::AbstractArray{<:Real}; kwargs...)
    return only(f(reshape(x, size(x, 1), size(x, 2), 1); kwargs...))
end

function _per_parameter(f, x::AbstractArray{<:AbstractArray{<:Real}}; kwargs...)
    shape = size(first(x))
    draws = permutedims(reshape(stack(x), prod(shape), size(x, 1), size(x, 2)), (2, 3, 1))
    return reshape(f(draws; kwargs...), shape)
end

end
