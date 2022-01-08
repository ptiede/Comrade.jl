export PolarizedModel, coherencymatrix

abstract type AbstractPolarizedModel end

"""
    $(TYPEDEF)
Wrapped model for a polarized model. This uses the stokes representation of the image.
"""
struct PolarizedModel{I,Q,U,V} <: AbstractPolarizedModel
    mi::I
    mq::Q
    mu::U
    mv::V
end


"""
    $(SIGNATURES)
Computes the visibility in the stokes basis of the polarized model
"""
@inline function visibility(pimg::PolarizedModel, u, v)
    i = visibility(pimg.i, u, v)
    q = visibility(pimg.q, u, v)
    u = visibility(pimg.u, u, v)
    v = visibility(pimg.v, u, v)
    return StokesVector(i, q, u, v)
end

"""
    $(SIGNATURES)
Finds the polarized intensity map of the polarized model pmodel.
"""
function intensitymap!(pimg::PolarizedMap, pmodel::PolarizedModel)
    intensitymap!(pimg.I, pmodel.I)
    intensitymap!(pimg.Q, pmodel.Q)
    intensitymap!(pimg.U, pmodel.U)
    intensitymap!(pimg.V, pmodel.V)
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
    return 1/2*atan2(su, sq)
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
