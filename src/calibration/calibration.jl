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
    JonesModel(jones::JonesPairs, tcache::TransformCache)

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

function intensitymap(model::RIMEModel, dims::AbstractDims)
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
    instrument = model.instrument
    jp = instrument.jones
    coh = _coherency(vis, typeof(instrument.refbasis))
    return corrupt(coh, jp.m1, jp.m2)
end

function ComradeBase.amplitudes(model::VLBIModel, ac::ArrayConfiguration)
    amp = amplitudes(model.sky, ac)
    instrument = model.instrument
    jp = instrument.jones
    coh = _coherency(amp, typeof(instrument.refbasis))
    return corrupt(coh, jp.m1, jp.m2)
end

function visibilities_numeric(model::VLBIModel, u, v, time, freq)
    vis = visibilities_numeric(model.sky, u, v, time, freq)
    instrument = model.instrument
    jp = instrument.jones
    coh = _coherency(vis, typeof(instrument.refbasis))
    return corrupt(coh, jp.m1, jp.m2)
end



include(joinpath(@__DIR__, "designmatrix.jl"))
include(joinpath(@__DIR__, "jones.jl"))
include(joinpath(@__DIR__, "priors.jl"))
include(joinpath(@__DIR__, "caltable.jl"))
