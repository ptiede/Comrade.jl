export PolarizedModel, coherencymatrix

import ComradeBase: AbstractPolarizedModel, m̆, evpa, CoherencyMatrix, StokesVector

"""
    $(TYPEDEF)
Wrapped model for a polarized model. This uses the stokes representation of the image.
"""
struct PolarizedModel{TI,TQ,TU,TV} <: AbstractPolarizedModel
    I::TI
    Q::TQ
    U::TU
    V::TV
end

@inline visanalytic(::Type{PolarizedModel{I,Q,U,V}}) where {I,Q,U,V} = visanalytic(I)*visanalytic(Q)*visanalytic(U)*visanalytic(V)
@inline imanalytic(::Type{PolarizedModel{I,Q,U,V}}) where {I,Q,U,V} = imanalytic(I)*imanalytic(Q)*imanalytic(U)*imanalytic(V)

@inline function intensity_point(pmodel::PolarizedModel, u, v)
    I = intensity_point(pmodel.I, u, v)
    Q = intensity_point(pmodel.Q, u, v)
    U = intensity_point(pmodel.U, u, v)
    V = intensity_point(pmodel.V, u, v)
    return StokesVector(I,Q,U,V)
end

"""
    $(SIGNATURES)
Computes the visibility in the stokes basis of the polarized model
"""
@inline function visibility(pimg::PolarizedModel, u, v)
    si = visibility(pimg.I, u, v)
    sq = visibility(pimg.Q, u, v)
    su = visibility(pimg.U, u, v)
    sv = visibility(pimg.V, u, v)
    return StokesVector(si, sq, su, sv)
end

"""
    $(SIGNATURES)
Finds the polarized intensity map of the polarized model pmodel.
"""
function intensitymap!(pimg::IntensityMap{<:StokesVector}, pmodel::PolarizedModel)
    intensitymap!(stokes(pimg, :I), pmodel.I)
    intensitymap!(stokes(pimg, :Q), pmodel.Q)
    intensitymap!(stokes(pimg, :U), pmodel.U)
    intensitymap!(stokes(pimg, :V), pmodel.V)
end

function intensitymap(pmodel::PolarizedModel, fovx::Real, fovy::Real, nx::Int, ny::Int; pulse=DeltaPulse())
    imgI = intensitymap(pmodel.I, fovx, fovy, nx, ny)
    imgQ = intensitymap(pmodel.Q, fovx, fovy, nx, ny)
    imgU = intensitymap(pmodel.U, fovx, fovy, nx, ny)
    imgV = intensitymap(pmodel.V, fovx, fovy, nx, ny)
    pimg = StructArray{StokesVector{eltype(imgI)}}((imgI.im, imgQ.im, imgU.im, imgV.im))
    return IntensityMap(pimg, fovx, fovy, pulse)
end

"""
    $(SIGNATURES)
Computes the coherency matrix of the polarized model.
"""
@inline function coherencymatrix(pimg::PolarizedModel, u, v)
    si = visibility(pimg, u, v)
    return convert(CoherencyMatrix, si)
end

"""
    $(SIGNATURES)
electric vector position angle or EVPA of the polarized model
"""
@inline function evpa(pimg, u, v)
    sq = visibility(pimg.Q, u, v)
    su = visibility(pimg.U, u, v)
    return 1/2*angle(su/sq)
end


"""
    $(SIGNATURES)
Computes the fractional linear polarization in the visibility domain
```math
    \\breve{m} = \\frac{Q + iU}{I}
```
"""
@inline function m̆(pimg, u, v)
    Q = visibility(pimg.Q, u, v)
    U = visibility(pimg.U, u, v)
    I = visibility(pimg.I, u, v)
    return (Q+1im*U)/I
end
