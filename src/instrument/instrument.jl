using LinearAlgebra
using SparseArrays
import Distributions
using Statistics
using PrettyTables

abstract type AbstractVLBIModel <: AbstractModel end
abstract type AbstractRIMEModel end


basemodel(m::AbstractVLBIModel) = m.sky
flux(m::AbstractVLBIModel) = flux(basemodel(m))
radialextent(m::AbstractVLBIModel) = radialextent(basemodel(m))

function intensitymap(model::AbstractVLBIModel, dims::AbstractGrid)
    return intensitymap(basemodel(model), dims)
end

function intensity_point(model::AbstractVLBIModel, p)
    return intensity_point(basemodel(model), p)
end

_amplitudes(m::M, u, v, t, f) where {M<:AbstractVLBIModel} = abs.(_visibilities(visanalytic(M), m, u, v, t, f))

"""
    VLBIModel(skymodel, instrumentmodel)

Constructs a `VLBIModel` from a `jones` pairs that describe the intrument model
and the `model` which describes the on-sky polarized visibilities. The third argument
can either be the `tcache` that converts from the model coherency basis to the instrumental
basis, or just the `refbasis` that will be used when constructing the model coherency matrices.
"""
struct VLBIModel{J<:AbstractRIMEModel, PI, M<:AbstractModel, PS} <: AbstractVLBIModel
    """
    The instrument model for the telescope. This is usually a sparse matrix that multiplies
    the visibilties.
    """
    instrument::J
    instrumentparams::PI
    """
    Base model that will be used to compute the uncorrupted Stokes visibilities.
    """
    sky::M
    skyparams::PS
end

visanalytic(::Type{<:VLBIModel{J,M}}) where {J,M} = visanalytic(M)
imanalytic(::Type{<:VLBIModel{J,M}}) where {J,M} = imanalytic(M)
ispolarized(::Type{<:VLBIModel{J,M}}) where {J,M} = ispolarized(M)


function visibilities_analytic(model::VLBIModel, u, v, time, freq)
    skym = model.sky(model.skyparams)
    vis = visibilities_analytic(skym, u, v, time, freq)
    return apply_instrument(vis, model.instrument, model.instrumentparams)
end

function visibilities_numeric(model::VLBIModel, u, v, time, freq)
    skym = model.sky
    vis = visibilities_numeric(skym, u, v, time, freq)
    return apply_instrument(vis, model.instrument)
end

struct IntegrationTime{T, I}
    mjd::T
    scanid::I
    t0::T
    dt::T
end

Base.in(t::Number, ts::IntegrationTime) = (ts.t0 - ts.dt/2) ≤ t < (ts.t0 + ts.dt/2)
scanid(ts::IntegrationTime) = ts.scanid
mjd(ts::IntegrationTime) = ts.mjd


struct FrequencyChannel{T, I<:Integer}
    central::T
    bandwidth::T
    channel::I
end

Base.in(f::Number, fs::FrequencyChannel) = (fs.central-fs.bandwidth/2) ≤ f < (fs.central+fs.bandwidth/2)
channel(fs::FrequencyChannel) = fs.channel


include("site_array.jl")
include("feedrotations.jl")
include("jonesmatrices.jl")
include("priors/priors.jl")
include("instrument_transforms.jl")
include("rime.jl")
include("utility.jl")
include("caltable.jl")
