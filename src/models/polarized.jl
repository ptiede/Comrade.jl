export PolarizedModel, coherencymatrix

import ComradeBase: AbstractPolarizedModel, m̆, evpa, CoherencyMatrix, StokesParams

"""
    $(TYPEDEF)

Wrapped model for a polarized model. This uses the stokes representation of the image.

# Fields
$(FIELDS)
"""
struct PolarizedModel{TI,TQ,TU,TV} <: AbstractPolarizedModel
    """
    Stokes I model
    """
    I::TI
    """
    Stokes Q Model
    """
    Q::TQ
    """
    Stokes U Model
    """
    U::TU
    """
    Stokes V Model
    """
    V::TV
end

function Base.show(io::IO, model::PolarizedModel)
    println(io, "PolarizedModel")
    println(io, "\tI: $(summary(model.I))")
    println(io, "\tQ: $(summary(model.Q))")
    println(io, "\tU: $(summary(model.U))")
    print(io, "\tV: $(summary(model.V))")
end

Base.@constprop :aggressive @inline visanalytic(::Type{PolarizedModel{I,Q,U,V}}) where {I,Q,U,V} = visanalytic(I)*visanalytic(Q)*visanalytic(U)*visanalytic(V)
Base.@constprop :aggressive @inline imanalytic(::Type{PolarizedModel{I,Q,U,V}}) where {I,Q,U,V} = imanalytic(I)*imanalytic(Q)*imanalytic(U)*imanalytic(V)

@inline function intensity_point(pmodel::PolarizedModel, u, v)
    I = intensity_point(pmodel.I, u, v)
    Q = intensity_point(pmodel.Q, u, v)
    U = intensity_point(pmodel.U, u, v)
    V = intensity_point(pmodel.V, u, v)
    return StokesParams(I,Q,U,V)
end

"""
    visibility(pimg::PolarizedModel, u, v)

Computes the visibility in the stokes basis of the polarized model
"""
@inline function visibility(pimg::PolarizedModel, u, v)
    si = visibility(pimg.I, u, v)
    sq = visibility(pimg.Q, u, v)
    su = visibility(pimg.U, u, v)
    sv = visibility(pimg.V, u, v)
    return StokesParams(si, sq, su, sv)
end

function visibilities(pimg::PolarizedModel, u, v)
    si = visibilities(pimg.I, u, v)
    sq = visibilities(pimg.Q, u, v)
    su = visibilities(pimg.U, u, v)
    sv = visibilities(pimg.V, u, v)
    return StructArray{StokesParams{eltype(si)}}((si, sq, su, sv))
end

function intensitymap!(pimg::StokesIntensityMap, pmodel::PolarizedModel)
    intensitymap!(stokes(pimg, :I), pmodel.I)
    intensitymap!(stokes(pimg, :Q), pmodel.Q)
    intensitymap!(stokes(pimg, :U), pmodel.U)
    intensitymap!(stokes(pimg, :V), pmodel.V)
    return pimg
end

function intensitymap(pmodel::PolarizedModel, dims)
    imgI = intensitymap(pmodel.I, dims)
    imgQ = intensitymap(pmodel.Q, dims)
    imgU = intensitymap(pmodel.U, dims)
    imgV = intensitymap(pmodel.V, dims)
    pimg = StructArray{StokesParams{eltype(imgI)}}((imgI.img, imgQ.img, imgU.img, imgV.img))
    return IntensityMap(pimg, fov, phasecenter, pulse)
end

# """
#     $(SIGNATURES)

# Computes the coherency matrix of the polarized model `pimg` at `u` and `v`
# """
# @inline function coherencymatrix(pimg::PolarizedModel, u, v)
#     si = visibility(pimg, u, v)
#     return convert(CoherencyMatrix, si)
# end

# """
#     evpa(pimg::AbstractPolarizedModel, u, v)

# electric vector position angle or EVPA of the polarized model `pimg` at `u` and `v`
# """
# @inline function evpa(pimg::AbstractPolarizedModel, u, v)
#     sq = visibility(pimg.Q, u, v)
#     su = visibility(pimg.U, u, v)
#     return 1/2*angle(su/sq)
# end


# """
#     m̆(pimg::AbstractPolarizedModel, u, v)

# Computes the fractional linear polarization in the visibility domain

#     m̆ = (Q + iU)/I

# To create the symbol type `m\\breve` in the REPL or use the
# [`mbreve`](@ref) function.
# """
# @inline function m̆(pimg::AbstractPolarizedModel, u, v)
#     Q = visibility(pimg.Q, u, v)
#     U = visibility(pimg.U, u, v)
#     I = visibility(pimg.I, u, v)
#     return (Q+1im*U)/I
# end

# """
#     $(SIGNATURES)

# Explicit m̆ function used for convenience.
# """
# mbreve(pimg, u, v) = m̆(pimg, u, v)
