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
    $(SIGNATURES)
Computes the complex visibility of model `m` at u,v positions `u,v`

If you want to compute the visibilities at a large number of positions
consider using the `visibilities` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function visibility(mimg::M, u, v) where {M}
    #first we split based on whether the model is primitive
    _visibility(isprimitive(M), mimg, u, v)
end

@inline function visibility(mimg, uv::ArrayBaselineDatum)
    return visibility(mimg, uv.u, uv.v)
end

"""
    $(SIGNATURES)
Computes the visibility amplitude of model `m` at u,v positions `u,v`

If you want to compute the amplitudes at a large number of positions
consider using the `amplitudes` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function amplitude(model, u, v)
    return abs(visibility(model, u, v))
end

"""
    $(SIGNATURES)
Computes the complex bispectrum of model `m` at the uv-triangle
u1,v1 -> u2,v2 -> u3,v3

If you want to compute the bispectrum over a number of triangles
consider using the `bispectra` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function bispectrum(model, u1, v1, u2, v2, u3, v3)
    return visibility(model, u1, v1)*visibility(model, u2, v2)*visibility(model, u3, v3)
end

"""
    $(SIGNATURES)
Computes the closure phase of model `m` at the uv-triangle
u1,v1 -> u2,v2 -> u3,v3

If you want to compute closure phases over a number of triangles
consider using the `closure_phases` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function closure_phase(model, u1, v1, u2, v2, u3, v3)
    return angle(bispectrum(model, u1, v1, u2, v2, u3, v3))
end

"""
    $(SIGNATURES)
Computes the log-closure amplitude of model `m` at the uv-quadrangle
u1,v1 -> u2,v2 -> u3,v3 -> u4,v3 using the formula

```math
C = \\log\\left|\\frac{V(u1,v1)V(u2,v2)}{V(u3,v3)V(u4,v4)}\\right|
```

If you want to compute log closure amplitudes over a number of triangles
consider using the `logclosure_amplitudes` function which uses MappedArrays to
compute a lazy, no allocation vector.
"""
@inline function logclosure_amplitude(model, u1, v1, u2, v2, u3, v3, u4, v4)
    a1 = amplitude(model, u1, v1)
    a2 = amplitude(model, u2, v2)
    a3 = amplitude(model, u3, v3)
    a4 = amplitude(model, u4, v4)

    return log(a1*a2/(a3*a4))
end


#=
    Welcome to the trait jungle. What is below is
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

@inline function visibilities(m::M, u::AbstractArray, v::AbstractArray) where {M}
    _visibilities(m, u, v)
    #visibilities(visanalytic(M), m, u, v)
end

function visibilities(m, ac::ArrayConfiguration)
    u, v = getuv(ac)
    return visibilities(m, u, v)
end

#isibilities(::IsAnalytic, m, u, v) = _visibilities(m, u, v)

#function create_mimg(m::AbstractModel, u, v)
#    fovx = fovy = radialextent(m)*2
#    ps = max(maximum(u), maximum(v))
#    nx, ny = abs(Int(ceil(10*fovx*ps))), abs(Int(ceil(10*fovy*ps)))
#    println(nx)
#   img = intensitymap(m, fovx, fovy, nx, ny)
#    return modelimage(m, img)
#end

#function create_mimg(m::AbstractModelImage, u, v)
#    return m
#end


#function visibilities(::NotAnalytic, m, u, v)
#   mimg = create_mimg(m, u, v)
#    println("here")
#    return _visibilities(mimg, u, v)
#end


"""
    $(SIGNATURES)
Computes the visibilities of the model `m` at the u,v positions `u`, `v`.

Note this is done lazily so the visibility is only computed when accessed.
"""
function _visibilities(m, u::AbstractArray, v::AbstractArray)
    f(x,y) = visibility(m, x, y)
    return mappedarray(f, u, v)
end



"""
    $(SIGNATURES)
Computes the amplitudes of the model `m` at the u,v positions `u`, `v`.

Note this is done lazily so the visibility is only computed when accessed.
"""
function amplitudes(m, u::AbstractArray, v::AbstractArray)
    f(x,y) = amplitude(m, x, y)
    return mappedarray(f, u, v)
end

"""
    $(SIGNATURES)
Computes the bispectra of the model `m` at the
triangles (u1,v1), (u2,v2), (u3,v3).

Note this is done lazily so the bispectra is only computed when accessed.
"""
function bispectra(m,
                    u1::AbstractArray, v1::AbstractArray,
                    u2::AbstractArray, v2::AbstractArray,
                    u3::AbstractArray, v3::AbstractArray
                   )
    f(x1,y1,x2,y2,x3,y3) = bispectra(m, x1, y1, x2, y2, x3, y3)
    return mappedarray(f, u1, v1, u2, v2, u3, v3)
end

"""
    $(SIGNATURES)
Computes the closure phases of the model `m` at the
triangles (u1,v1), (u2,v2), (u3,v3).

Note this is done lazily so the closure_phases is only computed when accessed.
"""
function closure_phases(m,
                        u1::AbstractArray, v1::AbstractArray,
                        u2::AbstractArray, v2::AbstractArray,
                        u3::AbstractArray, v3::AbstractArray
                       )
    f(x1,y1,x2,y2,x3,y3) = closure_phase(m, x1, y1, x2, y2, x3, y3)
    return mappedarray(f, u1, v1, u2, v2, u3, v3)
end

"""
    $(SIGNATURES)
Computes the log closure amplitudes of the model `m` at the
quadrangles (u1,v1), (u2,v2), (u3,v3), (u4, v4).

Note this is done lazily so the log closure amplitude is only computed when accessed.
"""
function logclosure_amplitudes(m,
                               u1::AbstractArray, v1::AbstractArray,
                               u2::AbstractArray, v2::AbstractArray,
                               u3::AbstractArray, v3::AbstractArray,
                               u4::AbstractArray, v4::AbstractArray
                              )
    f(x1,y1,x2,y2,x3,y3,x4,y4) = logclosure_amplitude(m, x1, y1, x2, y2, x3, y3, x4, y4)
    return mappedarray(f, u1, v1, u2, v2, u3, v3, u4, v4)
end



function intensitymap!(im::IntensityMap, m::M) where {M}
    return intensitymap!(imanalytic(M), im, m)
end


function intensitymap!(::NotAnalytic, img::IntensityMap, m)
    ny, nx = size(img)
    vis = ifftshift(phasedecenter!(fouriermap(m, img), img))
    ifft!(vis)
    for I in CartesianIndices(img)
        img[I] = real(vis[I])/(nx*ny)
    end
end


function intensitymap(::NotAnalytic, m, fovx::Real, fovy::Real, nx::Int, ny::Int; pulse=DeltaPulse())
    img = IntensityMap(zeros(ny, nx), fovx, fovy, pulse)
    vis = ifftshift(phasedecenter!(fouriermap(m, fovx, fovy, nx, ny), fovx, fovy, nx, ny))
    ifft!(vis)
    for I in CartesianIndices(img)
        img[I] = real(vis[I])/(nx*ny)
    end
    return img
end


include(joinpath(@__DIR__, "fft_alg.jl"))
include(joinpath(@__DIR__, "modelimage.jl"))
include(joinpath(@__DIR__, "modifiers.jl"))
include(joinpath(@__DIR__, "combinators.jl"))
include(joinpath(@__DIR__, "geometric_models.jl"))
include(joinpath(@__DIR__, "radio_image_models.jl"))
include(joinpath(@__DIR__, "polarized.jl"))
