export PolarizedModel, coherencymatrix, PoincareSphere2Map

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

@inline function intensity_point(pmodel::PolarizedModel, p)
    I = intensity_point(stokes(pmodel, :I), p)
    Q = intensity_point(stokes(pmodel, :Q), p)
    U = intensity_point(stokes(pmodel, :U), p)
    V = intensity_point(stokes(pmodel, :V), p)
    return StokesParams(I,Q,U,V)
end

ComradeBase.stokes(pmodel::PolarizedModel, symbol) = getproperty(pmodel, symbol)


"""
    visibility(pimg::PolarizedModel, p)

Computes the visibility in the stokes basis of the polarized model
"""
@inline function visibility(pimg::PolarizedModel, p)
    si = visibility(stokes(pimg, :I), p)
    sq = visibility(stokes(pimg, :Q), p)
    su = visibility(stokes(pimg, :U), p)
    sv = visibility(stokes(pimg, :V), p)
    return StokesParams(si, sq, su, sv)
end

function visibilities(pimg::PolarizedModel, p)
    si = visibilities(stokes(pimg, :I), p)
    sq = visibilities(stokes(pimg, :Q), p)
    su = visibilities(stokes(pimg, :U), p)
    sv = visibilities(stokes(pimg, :V), p)
    return StructArray{StokesParams{eltype(si)}}((si, sq, su, sv))
end

function intensitymap!(pimg::StokesIntensityMap, pmodel::PolarizedModel)
    intensitymap!(stokes(pimg, :I), pmodel.I)
    intensitymap!(stokes(pimg, :Q), pmodel.Q)
    intensitymap!(stokes(pimg, :U), pmodel.U)
    intensitymap!(stokes(pimg, :V), pmodel.V)
    return pimg
end

function intensitymap(pmodel::PolarizedModel, dims::Union{AbstractDims, NamedTuple})
    imgI = intensitymap(stokes(pmodel, :I), dims)
    imgQ = intensitymap(stokes(pmodel, :Q), dims)
    imgU = intensitymap(stokes(pmodel, :U), dims)
    imgV = intensitymap(stokes(pmodel, :V), dims)
    return StokesIntensityMap(imgI, imgQ, imgU, imgV)
end


"""
    PoincareSphere2Map(I, p, X, grid)
    PoincareSphere2Map(I::IntensityMap, p, X)

Constructs an polarized intensity map model using the Poincare parameterization.
The arguments are:
  - `I` is a grid of fluxes for each pixel.
  - `p` is a grid of numbers between 0, 1 and the represent the total fractional polarization
  - `X` is a grid, where each element is 3 numbers that represents the point on the Poincare sphere
    that is, X[1,1] is a NTuple{3} such that `||X[1,1]|| == 1`.
  - `grid` is the dimensional grid that gives the pixels locations of the intensity map.

!!! note
    If `I` is an `IntensityMap` then grid is not required since the same grid that was use
    for `I` will be used to construct the polarized intensity map

!!! warning
    The return type for this function is a polarized image object, however what we return
    is not considered to be part of the stable API so it may change suddenly.
"""
function PoincareSphere2Map(I, p, X, grid)
    pimgI = I.*p
    stokesI = IntensityMap(I, grid)
    stokesQ = IntensityMap(pimgI .* X[1], grid)
    stokesU = IntensityMap(pimgI .* X[2], grid)
    stokesV = IntensityMap(pimgI .* X[3], grid)
    return StokesIntensityMap(stokesI, stokesQ, stokesU, stokesV)
end
PoincareSphere2Map(I::IntensityMap, p, X) = PoincareSphere2Map(baseimage(I), p, X, ComradeBase.axiskeys(I))


"""
    evpa(pimg::AbstractPolarizedModel, p)

electric vector position angle or EVPA of the polarized model `pimg` at `u` and `v`
"""
@inline function evpa(pimg::AbstractPolarizedModel, p)
    sq = visibility(stokes(pimg :Q), p)
    su = visibility(stokes(pimg :U), p)
    return 1/2*angle(su/sq)
end


"""
    m̆(pimg::AbstractPolarizedModel, p)
    mbreve(pimg::AbstractPolarizedModel, p)

Computes the fractional linear polarization in the visibility domain

    m̆ = (Q + iU)/I

To create the symbol type `m\\breve` in the REPL or use the
[`mbreve`](@ref) function.
"""
@inline function m̆(pimg::AbstractPolarizedModel, p)
    Q = visibility(stokes(pimg, :Q), p)
    U = visibility(stokes(pimg, :U), p)
    I = visibility(stokes(pimg, :I), p)
    return (Q+1im*U)/I
end

"""
    $(SIGNATURES)

Explicit m̆ function used for convenience.
"""
mbreve(pimg, p) = m̆(pimg, p)
