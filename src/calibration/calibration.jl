using LinearAlgebra
using SparseArrays
import Distributions
using Statistics
using PrettyTables

export JonesModel

struct JonesModel{J, B}
    jones::J
    refbasis::B
end

"""
    JonesModel(jones::JonesPairs, refbasis = CirBasis())
    JonesModel(jones::JonesPairs, tcache::ResponseCache)

Constructs the intrument corruption model using pairs of jones matrices `jones` and a
reference basis
"""
JonesModel(jones::J) where {J} = JonesModel{J, CirBasis}(jones, CirBasis())



"""
    $(TYPEDEF)

Abstract type that encompasses all RIME style corruptions.
"""
abstract type RIMEModel <: AbstractModel end

basemodel(m::RIMEModel) = m.sky
flux(m::RIMEModel) = flux(basemodel(m))
radialextent(m::RIMEModel) = radialextent(basemodel(m))

function intensitymap(model::RIMEModel, dims::AbstractGrid)
    return intensitymap(basemodel(model), dims)
end

function intensity_point(model::RIMEModel, p)
    return intensity_point(basemodel(model), p)
end

_amplitudes(m::M, u, v, t, f) where {M<:RIMEModel} = abs.(_visibilities(visanalytic(M), m, u, v, t, f))

"""
    VLBIModel(skymodel, instrumentmodel)

Constructs a `VLBIModel` from a `jones` pairs that describe the intrument model
and the `model` which describes the on-sky polarized visibilities. The third argument
can either be the `tcache` that converts from the model coherency basis to the instrumental
basis, or just the `refbasis` that will be used when constructing the model coherency matrices.
"""
struct VLBIModel{J<:JonesModel, M<:AbstractModel} <: RIMEModel
    """
    The instrument model for the telescope. This is usually a sparse matrix that multiplies
    the visibilties.
    """
    instrument::J
    """
    Base model that will be used to compute the uncorrupted Stokes visibilities.
    """
    sky::M
end

visanalytic(::Type{<:VLBIModel{J,M}}) where {J,M} = visanalytic(M)
imanalytic(::Type{<:VLBIModel{J,M}}) where {J,M} = imanalytic(M)
ispolarized(::Type{<:VLBIModel{J,M}}) where {J,M} = ispolarized(M)


function visibilities_analytic(model::VLBIModel, u, v, time, freq)
    vis = visibilities_analytic(model.sky, u, v, time, freq)
    return apply_instrument(vis, model.instrument)
end

function visibilities_numeric(model::VLBIModel, u, v, time, freq)
    skym = model.sky
    vis = visibilities_numeric(skym, u, v, time, freq)
    return apply_instrument(vis, model.instrument)
end

function ComradeBase.amplitudes(model::VLBIModel, ac::ArrayConfiguration)
    amp = amplitudes(model.sky, ac)
    apply_instrument(amp, model.instrument)
end

function apply_instrument(vis, instrument)
    jp = instrument.jones
    vout = _apply_instrument(vis, jp.m1, jp.m2, instrument.refbasis)
    return vout
end

intout(vis::AbstractArray{<:StokesParams{T}}, j1, j2) where {T<:Real} = similar(vis, SMatrix{2,2, Complex{T}, 4})
intout(vis::AbstractArray{T}, j1, j2) where {T<:Real} = similar(vis, Complex{T})
intout(vis::AbstractArray{<:CoherencyMatrix{A,B,T}}, j1, j2) where {A,B,T<:Real} = similar(vis, SMatrix{2,2, Complex{T}, 4})

intout(vis::AbstractArray{<:StokesParams{T}}, j1, j2) where {T<:Complex} = similar(vis, SMatrix{2,2, T, 4})
intout(vis::AbstractArray{T}, j1, j2) where {T<:Complex} = similar(vis, T)
intout(vis::AbstractArray{<:CoherencyMatrix{A,B,T}}, j1, j2) where {A,B,T<:Complex} = similar(vis, SMatrix{2,2, T, 4})


function _apply_instrument(vis, j1, j2, refbasis)
    vout = intout(vis, j1, j2)
    _apply_instrument!(vout, vis, j1, j2, refbasis)
    return vout
end

function _apply_instrument!(vout, vis, jp1, jp2, refbasis)
    map!(vout, vis, jp1, jp2) do vi, j1, j2
        _apply_jones(vi, j1, j2, refbasis)
    end
    return nothing
end

_apply_jones(v::Number, j1, j2, ::B) where {B} = j1*v*j2'
_apply_jones(v::CoherencyMatrix, j1, j2, ::B) where {B} = j1*CoherencyMatrix{B,B}(v)*j2'
_apply_jones(v::StokesParams, j1, j2, ::B) where {B} = j1*CoherencyMatrix{B,B}(v)*j2'


function ChainRulesCore.rrule(::typeof(_apply_instrument), vis, j1, j2, refbasis)
    out = _apply_instrument(vis, j1, j2, refbasis)
    pvis = ProjectTo(vis)
    pj1 = ProjectTo(j1)
    pj2 = ProjectTo(j2)
    function _apply_instrument_pullback(Δ)
        Δvout = zero(out)
        Δvout .= unthunk(Δ)
        djp1 = zero(j1)
        djp2 = zero(j2)
        dvis = zero(vis)

        vout = similar(out)
        autodiff(Reverse, _apply_instrument!, Const,
                        Duplicated(vout, Δvout),
                        Duplicated(vis, dvis),
                        Duplicated(j1, djp1),
                        Duplicated(j2, djp2),
                        Const(refbasis)
                )

        return NoTangent(), pvis(dvis), pj1(djp1), pj2(djp2), NoTangent()
    end
    return out, _apply_instrument_pullback
end






include(joinpath(@__DIR__, "designmatrix.jl"))
include(joinpath(@__DIR__, "jones.jl"))
include(joinpath(@__DIR__, "priors.jl"))
include(joinpath(@__DIR__, "caltable.jl"))
