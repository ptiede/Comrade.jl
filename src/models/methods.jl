@inline function ComradeBase.visibilities(m::M, p::ArrayConfiguration) where {M <: AbstractModel}
    return _visibilities(visanalytic(M), m, p.data.U, p.data.V, p.data.T, p.data.F)
end

@inline function ComradeBase.visibilities(m::M, p::ClosureConfig) where {M <: AbstractModel}
    return visibilities(m, arrayconfig(p.ac))
end


@inline function ComradeBase.amplitudes(m::AbstractModel, p::ArrayConfiguration)
    return amplitudes(m, (U = p.data.U, V = p.data.V, T=p.data.T, F=p.data.F))
end


"""
    closure_phases(m::AbstractModel, ac::ClosureConfig)

Computes the closure phases of the model `m` using the array configuration `ac`.

# Notes
This is faster than the `closure_phases(m, u1, v1, ...)` method since it only
computes as many visibilities as required thanks to the closure design matrix formalism
from Blackburn et al.[^1]

[^1]: Blackburn L., et al "Closure Statistics in Interferometric Data" ApJ 2020
"""
function ComradeBase.closure_phases(m::AbstractModel, ac::ClosureConfig)
    vis = visibilities(m, arrayconfig(ac.ac))
    return ac.designmat*angle.(vis)
end

"""
    logclosure_amplitudes(m::AbstractModel, ac::ClosureConfig)

Computes the log closure amplitudes of the model `m` using the array configuration `ac`.

# Notes
This is faster than the `logclosure_amplitudes(m, u1, v1, ...)` method since it only
computes as many visibilities as required thanks to the closure design matrix formalism
from Blackburn et al.[^1]

[^1]: Blackburn L., et al "Closure Statistics in Interferometric Data" ApJ 2020
"""
function ComradeBase.logclosure_amplitudes(m::AbstractModel, ac::ClosureConfig)
    vis = visibilities(m, arrayconfig(ac.ac))
    return ac.designmat*log.(abs.(vis))
end
