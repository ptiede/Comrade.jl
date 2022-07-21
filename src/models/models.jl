import ComradeBase: AbstractModel, IsPrimitive, NotPrimitive, IsAnalytic, NotAnalytic,
                    visanalytic, imanalytic, isprimitive
import ComradeBase: visibility_point,
                    intensitymap, intensitymap!, intensity_point,
                    flux

export visibility, amplitude, closure_phase, logclosure_amplitude, bispectrum,
       visibilities, amplitudes, closure_phases, logclosure_amplitudes, bispectra,
       flux, intensitymap, intensitymap!, PolarizedModel


abstract type AbstractModelImage{M} <: ComradeBase.AbstractModel end


"""
    visibility(mimg, u, v, args...)

Computes the complex visibility of model `m` at u,v positions `u,v` and
other `args...` (e.g., frequency, time, etc.)

# Notes
If you want to compute the visibilities at a large number of positions
consider using the [`visibilities`](@ref visibilities).
"""
@inline function visibility(mimg::M, u, v, args...) where {M}
    #first we split based on whether the model is primitive
    _visibility(isprimitive(M), mimg, u, v, args...)
end


"""
    visibility(mimg, uv::ArrayBaselineDatum)

Computes the complex visibility of a model `m` at the `uv` array baseline datum.

# Notes
If you want to compute the visibilities at a large number of positions
consider using the [`visibilities`](@ref visibilities).
"""
@inline function visibility(mimg, uv::ArrayBaselineDatum)
    return visibility(mimg, uv.u, uv.v, )
end

"""
    amplitude(model, args...)

Computes the visibility amplitude of model `m` at u,v positions `u,v`

If you want to compute the amplitudes at a large number of positions
consider using the `amplitudes` function.
"""
@inline function amplitude(model, args...)
    return abs(visibility(model, args...))
end

"""
    bispectrum(model, u1, v1, u2, v2, u3, v3)

Computes the complex bispectrum of model `m` at the uv-triangle
u1,v1 -> u2,v2 -> u3,v3

If you want to compute the bispectrum over a number of triangles
consider using the `bispectra` function.
"""
@inline function bispectrum(model, u1, v1, u2, v2, u3, v3)
    return visibility(model, u1, v1)*visibility(model, u2, v2)*visibility(model, u3, v3)
end

"""
    closure_phase(model. u1, v1, u2, v2, u3, v3)

Computes the closure phase of model `m` at the uv-triangle
u1,v1 -> u2,v2 -> u3,v3

If you want to compute closure phases over a number of triangles
consider using the `closure_phases` function.
"""
@inline function closure_phase(model, u1, v1, u2, v2, u3, v3)
    return angle(bispectrum(model, u1, v1, u2, v2, u3, v3))
end

"""
    logclosure_amplitude(model, u1, v1, u2, v2, u3, v3, u4, v4)

Computes the log-closure amplitude of model `m` at the uv-quadrangle
u1,v1 -> u2,v2 -> u3,v3 -> u4,v4 using the formula

```math
C = \\log\\left|\\frac{V(u1,v1)V(u2,v2)}{V(u3,v3)V(u4,v4)}\\right|
```

If you want to compute log closure amplitudes over a number of triangles
consider using the `logclosure_amplitudes` function.
"""
@inline function logclosure_amplitude(model, u1, v1, u2, v2, u3, v3, u4, v4)
    a1 = amplitude(model, u1, v1)
    a2 = amplitude(model, u2, v2)
    a3 = amplitude(model, u3, v3)
    a4 = amplitude(model, u4, v4)

    return log(a1*a2/(a3*a4))
end


#=
    Welcome to the trait jungle. Below is
    how we specify how to evaluate the model
=#
@inline function _visibility(::NotPrimitive, m, u, v)
    return visibility_point(m, u, v)
end

@inline function _visibility(::IsPrimitive, m::M, u, v) where {M}
    _visibility_primitive(visanalytic(M), m, u, v)
end


@inline function _visibility_primitive(::IsAnalytic, mimg, u, v)
    return visibility_point(mimg, u, v)
end

@inline function _visibility_primitive(::NotAnalytic, mimg, u, v)
    return mimg.cache.sitp(u, v)
end

# @inline function visibilities(m::M, u::AbstractArray, v::AbstractArray) where {M}
#     _visibilities(m, u, v)
# end

"""
    visibilities(m, ac::ArrayConfiguration)

Computes the visibilities of the model `m` using the array configuration `ac`.
If you want to compute a single visibility you should call [`visibility`](@ref visibility).
"""
@inline function visibilities(m, ac::ArrayConfiguration)
    u, v = getuv(ac)
    return visibilities(m, u, v)
end

"""
    visibilities(m, u::AbstractArray, v::AbstractArray, args...)

Computes the visibilities of the model `m` at `u` `v` and `args...`.
If you want to compute a single visibility you should call [`visibility`](@ref visibility).
"""
@inline function visibilities(m, u::AbstractArray, v::AbstractArray, args...)
    return _visibilities(m, u, v, args...)
end



# Internal function required for dispatch. This is a fallback method if
# visibilities doesn't have a direct implementation.
@inline function _visibilities(m, u::AbstractArray, v::AbstractArray, args...)
    _visibilities_fallback(m, u, v, args...)
end

function _visibilities_fallback(m, u::AbstractArray, v::AbstractArray, args...)
    #vm(u,v) = visibility(m, u, v)
    #return map(vm, u, v)
    return visibility.(Ref(m), u, v, args...)
end




"""
    amplitudes(m::AbstractModel, u::AbstractArray, v::AbstractArray)

Computes the amplitudes of the model `m` at the u,v positions `u`, `v`.
"""
function amplitudes(m, u::AbstractArray, v::AbstractArray)
    _amplitudes(m, u, v)
end

function amplitudes(m, ac::ArrayConfiguration)
    u, v = getuv(ac)
    return amplitudes(m, u, v)
end

function _amplitudes(m::S, u::AbstractArray, v::AbstractArray) where {S}
    _amplitudes(visanalytic(S), m, u, v)
end

function _amplitudes(::IsAnalytic, m, u::AbstractArray, v::AbstractArray)
    #f(x,y) = amplitude(m, x, y)
    return amplitude.(Ref(m), u, v)
end

function _amplitudes(::NotAnalytic, m, u::AbstractArray, v::AbstractArray)
    abs.(visibilities(m, u, v))
end


"""
    bispectra(m,
              u1::AbstractArray, v1::AbstractArray,
              u2::AbstractArray, v2::AbstractArray,
              u3::AbstractArray, v3::AbstractArray,
            )

Computes the bispectra of the model `m` at the
triangles (u1,v1), (u2,v2), (u3,v3).
"""
function bispectra(m,
                    u1::AbstractArray, v1::AbstractArray,
                    u2::AbstractArray, v2::AbstractArray,
                    u3::AbstractArray, v3::AbstractArray
                    )
    _bispectra(m, u1, v1, u2, v2, u3, v3)
end

# internal method used for trait dispatch
function _bispectra(m::M,
                    u1::AbstractArray, v1::AbstractArray,
                    u2::AbstractArray, v2::AbstractArray,
                    u3::AbstractArray, v3::AbstractArray
                    ) where {M}
    _bispectra(visanalytic(M), m, u1, v1, u2, v2, u3, v3)
end

# internal method used for trait dispatch for analytic visibilities
function _bispectra(::IsAnalytic, m,
                    u1::AbstractArray, v1::AbstractArray,
                    u2::AbstractArray, v2::AbstractArray,
                    u3::AbstractArray, v3::AbstractArray
                   )
    return bispectrum.(Ref(m), u1, v1, u2, v2, u3, v3)
end

# internal method used for trait dispatch for non-analytic visibilities
function _bispectra(::NotAnalytic, m,
                    u1::AbstractArray, v1::AbstractArray,
                    u2::AbstractArray, v2::AbstractArray,
                    u3::AbstractArray, v3::AbstractArray
                   )
    vis1 = _visibilities(m, u1, v1)
    vis2 = _visibilities(m, u2, v2)
    vis3 = _visibilities(m, u3, v3)
    return @. vis1*vis2*vis3
end

"""
    closure_phases(m,
                   u1::AbstractArray, v1::AbstractArray,
                   u2::AbstractArray, v2::AbstractArray,
                   u3::AbstractArray, v3::AbstractArray,
                   )

Computes the closure phases of the model `m` at the
triangles (u1,v1), (u2,v2), (u3,v3).
"""
@inline function closure_phases(m::AbstractModel,
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray
                        )
    _closure_phases(m, u1, v1, u2, v2, u3, v3)
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
@inline function _closure_phases(m::M,
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray
                       ) where {M<:AbstractModel}
    _closure_phases(visanalytic(M), m, u1, v1, u2, v2, u3, v3)
end

# internal method used for trait dispatch for analytic visibilities
@inline function _closure_phases(::IsAnalytic, m,
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray
                       )
    return closure_phase.(Ref(m), u1, v1, u2, v2, u3, v3)
end

# internal method used for trait dispatch for non-analytic visibilities
function _closure_phases(::NotAnalytic, m,
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray
                       )
    return angle.(bispectra(m, u1, v1, u2, v2, u3, v3))
end

"""
    logclosure_amplitudes(m::AbstractModel,
                          u1::AbstractArray, v1::AbstractArray,
                          u2::AbstractArray, v2::AbstractArray,
                          u3::AbstractArray, v3::AbstractArray,
                          u4::AbstractArray, v4::AbstractArray,
                         )

Computes the log closure amplitudes of the model `m` at the
quadrangles (u1,v1), (u2,v2), (u3,v3), (u4, v4).
"""
function logclosure_amplitudes(m::AbstractModel,
                               u1::AbstractArray, v1::AbstractArray,
                               u2::AbstractArray, v2::AbstractArray,
                               u3::AbstractArray, v3::AbstractArray,
                               u4::AbstractArray, v4::AbstractArray
                              )
    _logclosure_amplitudes(m, u1, v1, u2, v2, u3, v3, u4, v4)
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
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray,
                        u4::AbstractArray, v4::AbstractArray
                       ) where {M<:AbstractModel}
    _logclosure_amplitudes(visanalytic(M), m, u1, v1, u2, v2, u3, v3, u4, v4)
end

# internal method used for trait dispatch for analytic visibilities
@inline function _logclosure_amplitudes(::IsAnalytic, m,
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray,
                        u4::AbstractArray, v4::AbstractArray)

    return logclosure_amplitude.(Ref(m), u1, v1, u2, v2, u3, v3, u4, v4)
end

# internal method used for trait dispatch for non-analytic visibilities
@inline function _logclosure_amplitudes(::NotAnalytic, m,
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray,
                        u4::AbstractArray, v4::AbstractArray
                       )
    amp1 = amplitudes(m, u1, v1)
    amp2 = amplitudes(m, u2, v2)
    amp3 = amplitudes(m, u3, v3)
    amp4 = amplitudes(m, u4, v4)
    return @. log(amp1*amp2*inv(amp3*amp4))
end


# internal method for computing an image of a non-analytic image model. The
# `executor` if for parallelization but is not used for this method.
function intensitymap!(::NotAnalytic, img::IntensityMap, m, executor=SequentialEx())
    ny, nx = size(img)
    fovx, fovy = fov(img)
    vis = fouriermap(m, fovx, fovy, nx, ny)
    vis = ifftshift(phasedecenter!(vis, fovx, fovy, nx, ny))
    ifft!(vis)
    for I in CartesianIndices(img)
        img[I] = real(vis[I])/(nx*ny)
    end
end

# internal method for computing an image of a non-analytic image model. The
# `executor` if for parallelization but is not used for this method.
function intensitymap(::NotAnalytic, m, fovx::Real, fovy::Real, nx::Int, ny::Int; pulse=DeltaPulse(), executor=SequentialEx())
    img = IntensityMap(zeros(ny, nx), fovx, fovy, pulse)
    vis = ifftshift(phasedecenter!(fouriermap(m, fovx, fovy, nx, ny), fovx, fovy, nx, ny))
    ifft!(vis)
    for I in CartesianIndices(img)
        img[I] = real(vis[I])/(nx*ny)
    end
    return img
end


include(joinpath(@__DIR__, "modelimage/modelimage.jl"))
include(joinpath(@__DIR__, "modifiers.jl"))
include(joinpath(@__DIR__, "combinators.jl"))
include(joinpath(@__DIR__, "geometric_models.jl"))
include(joinpath(@__DIR__, "polarized.jl"))
include(joinpath(@__DIR__, "test.jl"))
include(joinpath(@__DIR__, "misc.jl"))
