"""
    visibility(mimg, p)

Computes the complex visibility of model `m` at coordinates `p`. `p` corresponds to
the coordinates of the model. These need to have the properties `U`, `V` and sometimes
`Ti` for time and `Fr` for frequency.

# Notes
If you want to compute the visibilities at a large number of positions
consider using the [`visibilities`](@ref visibilities).
"""
@inline function visibility(mimg::M, p) where {M}
    #first we split based on whether the model is primitive
    _visibility(isprimitive(M), mimg, p.U, p.V, 0.0, 0.0)
end


"""
    amplitude(model, p)

Computes the visibility amplitude of model `m` at the coordinate `p`.
The coordinate `p`
is expected to have the properties `U`, `V`, and sometimes `Ti` and `Fr`.

If you want to compute the amplitudes at a large number of positions
consider using the `amplitudes` function.
"""
@inline function amplitude(model, p)
    return abs(visibility(model, p))
end

"""
    bispectrum(model, p1, p2, p3)

Computes the complex bispectrum of model `m` at the uv-triangle
p1 -> p2 -> p3

If you want to compute the bispectrum over a number of triangles
consider using the `bispectra` function.
"""
@inline function bispectrum(model, p1, p2, p3)
    return visibility(model, p1)*visibility(model, p2)*visibility(model, p3)
end

"""
    closure_phase(model, p1, p2, p3, p4)

Computes the closure phase of model `m` at the uv-triangle
u1,v1 -> u2,v2 -> u3,v3

If you want to compute closure phases over a number of triangles
consider using the `closure_phases` function.
"""
@inline function closure_phase(model, p1, p2, p3)
    return angle(bispectrum(model, p1, p2, p3))
end

"""
    logclosure_amplitude(model, p1, p2, p3, p4)

Computes the log-closure amplitude of model `m` at the uv-quadrangle
u1,v1 -> u2,v2 -> u3,v3 -> u4,v4 using the formula

```math
C = \\log\\left|\\frac{V(u1,v1)V(u2,v2)}{V(u3,v3)V(u4,v4)}\\right|
```

If you want to compute log closure amplitudes over a number of triangles
consider using the `logclosure_amplitudes` function.
"""
@inline function logclosure_amplitude(model, p1, p2, p3, p4)
    a1 = amplitude(model, p1)
    a2 = amplitude(model, p2)
    a3 = amplitude(model, p3)
    a4 = amplitude(model, p4)

    return log(a1*a2/(a3*a4))
end


#=
    Welcome to the trait jungle. Below is
    how we specify how to evaluate the model
=#
@inline function _visibility(::NotPrimitive, m, u, v, time, freq)
    return visibility_point(m, u, v, time, freq)
end

@inline function _visibility(::IsPrimitive, m::M, u, v, time, freq) where {M}
    _visibility_primitive(visanalytic(M), m, u, v, time, freq)
end


@inline function _visibility_primitive(::IsAnalytic, mimg, u, v, time, freq)
    return visibility_point(mimg, u, v, time, freq)
end

@inline function _visibility_primitive(::NotAnalytic, mimg, u, v, time, freq)
    return mimg.cache.sitp(u, v)
end

# @inline function visibilities(m::M, u::AbstractArray, v::AbstractArray) where {M}
#     _visibilities(m, u, v)
# end

"""
    visibilities(m, p)

Computes the visibilities of the model `m` using the coordinates `p`. The coordinates `p`
are expected to have the properties `U`, `V`, and sometimes `Ti` and `Fr`.
"""
@inline function visibilities(m::AbstractModel, p::NamedTuple)
    U, V, T, F = extract_pos(p)
    return _visibilities(m, U, V, T, F)
end

function extract_pos(p::NamedTuple)
    return p.U, p.V, p.T, p.F
end

function extract_pos(p::NamedTuple{(:U,:V)})
    return p.U, p.V, zero(eltype(p.U)), zero(eltype(p.V))
end

# @inline function visibilities(m::AbstractModel, p::NamedTuple{(:U, :V)})
#     return _visibilities(m, p.U, p.V, zero(eltype(p.U)), zero(eltype(p.U)))
# end


@inline function visibilities(m, p::ArrayConfiguration)
    return _visibilities(m, p.data.U, p.data.V, p.data.T, p.data.F)
end

# function ChainRulesCore.rrule(::typeof(getuv), p::ArrayConfiguration)
#     sa = getuv(p)
#     function _get_uv_pullback(Δ)
#         println(typeof(Δ))
#         return Tangent{typeof(p)}(data = Tangent{typeof(p.data)}(U=Δ.U, V=Δ.V))
#     end
#     return sa, _get_uv_pullback
# end




# Internal function required for dispatch. This is a fallback method if
# visibilities doesn't have a direct implementation.
@inline function _visibilities(m, u, v, time, freq)
    _visibilities_fallback(m, u, v, time, freq)
end

function _visibilities_fallback(m, u, v, time, freq)
    return visibility_point.(Ref(m), u, v, time, freq)
end





"""
    amplitudes(m::AbstractModel, u::AbstractArray, v::AbstractArray)

Computes the visibility amplitudes of the model `m` at the coordinates `p`.
The coordinates `p` are expected to have the properties `U`, `V`,
and sometimes `Ti` and `Fr`.
"""
function amplitudes(m, p::NamedTuple{(:U, :V, :T, :F)})
    _amplitudes(m, p.U, p.V, p.T, p.F)
end

function amplitudes(m, p::NamedTuple{(:U, :V)})
    _amplitudes(m, p.U, p.V, 0.0, 0.0)
end

amplitudes(m, p::ArrayConfiguration) = _amplitudes(m, p.data.U, p.data.V, p.data.T, p.data.F)


function _amplitudes(m::S, u, v, time, freq) where {S}
    _amplitudes(visanalytic(S), m, u, v, time, freq)
end

function _amplitudes(::IsAnalytic, m, u, v, time, freq)
    abs.(visibility_point.(Ref(m), u, v, time, freq))
end

function _amplitudes(::NotAnalytic, m, u, v, time, freq)
    abs.(_visibilities(m, u, v, time, freq))
end


"""
    bispectra(m, p1, p2, p3)

Computes the closure phases of the model `m` at the
triangles p1, p2, p3, where `pi` are coordinates.
"""
function bispectra(m,
                    p1,
                    p2,
                    p3,
                    )

    _bispectra(m, p1, p2, p3)
end

# internal method used for trait dispatch
function _bispectra(m::M,
                    p1,
                    p2,
                    p3
                    ) where {M}
    _bispectra(visanalytic(M), m, p1, p2, p3)
end

# internal method used for trait dispatch for analytic visibilities
function _bispectra(::IsAnalytic, m,
                    p1,
                    p2,
                    p3,
                   )
    return bispectrum.(Ref(m), StructArray(p1), StructArray(p2), StructArray(p3))
end

# internal method used for trait dispatch for non-analytic visibilities
function _bispectra(::NotAnalytic, m,
                    p1,p2,p3
                   )
    vis1 = visibilities(m, p1)
    vis2 = visibilities(m, p2)
    vis3 = visibilities(m, p3)
    return @. vis1*vis2*vis3
end

"""
    closure_phases(m,
                   p1::AbstractArray
                   p2::AbstractArray
                   p3::AbstractArray
                   )

Computes the closure phases of the model `m` at the
triangles p1, p2, p3, where `pi` are coordinates.
"""
@inline function closure_phases(m::AbstractModel,
                        p1,p2,p3
                        )
    _closure_phases(m, p1, p2, p3)
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
function closure_phases(m::AbstractModel, ac::ClosureConfig)
    vis = visibilities(m, arrayconfig(ac.ac))
    return ac.designmat*angle.(vis)
end

# internal method used for trait dispatch
@inline function _closure_phases(m::M, p1, p2, p3) where {M<:AbstractModel}
    _closure_phases(visanalytic(M), m, p1, p2, p3)
end

# internal method used for trait dispatch for analytic visibilities
@inline function _closure_phases(::IsAnalytic, m,
                        p1,
                        p2,
                        p3
                       )
    return closure_phase.(Ref(m), p1, p2, p3)
end

# internal method used for trait dispatch for non-analytic visibilities
function _closure_phases(::NotAnalytic, m, p1,p2, p3)
    return angle.(bispectra(m, p1, p2, p3))
end

"""
    logclosure_amplitudes(m::AbstractModel,
                          p1,
                          p2,
                          p3,
                          p4
                         )

Computes the log closure amplitudes of the model `m` at the
quadrangles p1, p2, p3, p4.
"""
function logclosure_amplitudes(m::AbstractModel,
                               p1,
                               p2,
                               p3,
                               p4
                              )
    _logclosure_amplitudes(m, p1, p2, p3, p4)
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
function logclosure_amplitudes(m::AbstractModel, ac::ClosureConfig)
    vis = visibilities(m, arrayconfig(ac.ac))
    return ac.designmat*log.(abs.(vis))
end

# internal method used for trait dispatch
@inline function _logclosure_amplitudes(m::M,
                        p1,
                        p2,
                        p3,
                        p4
                       ) where {M<:AbstractModel}
    _logclosure_amplitudes(visanalytic(M), m, p1, p2, p3, p4)
end

# internal method used for trait dispatch for analytic visibilities
@inline function _logclosure_amplitudes(::IsAnalytic, m, p1, p2, p3, p4)

    return logclosure_amplitude.(Ref(m), p1, p2, p3, p4)
end

# internal method used for trait dispatch for non-analytic visibilities
@inline function _logclosure_amplitudes(::NotAnalytic, m, p1, p2, p3, p4)
    amp1 = amplitudes(m, p1)
    amp2 = amplitudes(m, p2)
    amp3 = amplitudes(m, p3)
    amp4 = amplitudes(m, p4)
    return @. log(amp1*amp2*inv(amp3*amp4))
end


# internal method for computing an image of a non-analytic image model. The
# `executor` if for parallelization but is not used for this method.
function intensitymap!(::NotAnalytic, img::IntensityMap, m)
    # nx, ny = size(img)
    (;X, Y) = axisdims(img)
    vis = fouriermap(m, axisdims(img))
    vis = ifftshift(phasedecenter!(vis, X, Y))
    ifft!(vis)
    img .= real.(vis[I])
    return
end

function intensitymap(A::NotAnalytic, m, grid::AbstractDims)
    img = IntensityMap(zeros(map(length, dims(grid))), grid)
    intensitymap!(A, img, m)
    return img
end


# function intensitymap(::NotAnalytic, m, dims)
#     vis = ifftshift(ComradeBase.AxisKeys.keyless_unname(phasedecenter!(fouriermap(m, dims), dims.X, dims.Y)))
#     ifft!(vis)
#     return IntensityMap(real.(vis)./length(vis), dims)
# end
