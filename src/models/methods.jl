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
    _visibility(isprimitive(M), mimg, p)
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
@inline function _visibility(::NotPrimitive, m, p)
    return visibility_point(m, p)
end

@inline function _visibility(::IsPrimitive, m::M, p) where {M}
    _visibility_primitive(visanalytic(M), m, p)
end


@inline function _visibility_primitive(::IsAnalytic, mimg, p)
    return visibility_point(mimg, p)
end

@inline function _visibility_primitive(::NotAnalytic, mimg, p)
    return mimg.cache.sitp(p.U, p.V)
end

# @inline function visibilities(m::M, u::AbstractArray, v::AbstractArray) where {M}
#     _visibilities(m, u, v)
# end

"""
    visibilities(m, p)

Computes the visibilities of the model `m` using the coordinates `p`. The coordinates `p`
are expected to have the properties `U`, `V`, and sometimes `Ti` and `Fr`.
"""
@inline function visibilities(m, p)
    return _visibilities(m, p)
end

@inline function visibilities(m, p::ArrayConfiguration)
    return _visibilities(m, getuv(p))
end




# Internal function required for dispatch. This is a fallback method if
# visibilities doesn't have a direct implementation.
@inline function _visibilities(m, p)
    _visibilities_fallback(m, p)
end

function _visibilities_fallback(m, p::NamedTuple)
    return visibility_point.(Ref(m), NamedTuple{(:U, :V)}.(p.U, p.V))
end

function _visibilities_fallback(m, p::StructArray)
    return visibility_point.(Ref(m), NamedTuple{(:U, :V)}.(p.U, p.V))
end





"""
    amplitudes(m::AbstractModel, u::AbstractArray, v::AbstractArray)

Computes the visibility amplitudes of the model `m` at the coordinates `p`.
The coordinates `p` are expected to have the properties `U`, `V`,
and sometimes `Ti` and `Fr`.
"""
function amplitudes(m, p)
    _amplitudes(m, p)
end

function _amplitudes(m::S, p) where {S}
    _amplitudes(visanalytic(S), m, p)
end

function _amplitudes(::IsAnalytic, m, p)
    #f(x,y) = amplitude(m, x, y)
    return amplitude.(Ref(m), p)
end

function _amplitudes(::NotAnalytic, m, p)
    abs.(visibilities(m, p))
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
    return bispectrum.(Ref(m), p1, p2, p3)
end

# internal method used for trait dispatch for non-analytic visibilities
function _bispectra(::NotAnalytic, m,
                    p1,p2,p3
                   )
    vis1 = _visibilities(m, p1)
    vis2 = _visibilities(m, p2)
    vis3 = _visibilities(m, p3)
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
    vis = visibilities(m, ac.ac)
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
    vis = visibilities(m, ac.ac)
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
    ny, nx = size(img)
    vis = fouriermap(m, dims(img))
    vis = ifftshift(phasedecenter!(vis, xitr, yitr))
    ifft!(vis, dimnum(img, (X,Y)))
    img .= real.(vis[I])./(nx*ny)

end

function intensitymap(::NotAnalytic, m, dims)
    vis = ifftshift(phasedecenter!(fouriermap(m, dims), dims.X, dims.Y))
    ifft!(vis)
    return real.(vis)./length(vis)
end
