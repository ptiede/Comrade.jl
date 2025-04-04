using LinearAlgebra
using SparseArrays
import Distributions
using Statistics
using PrettyTables

struct IntegrationTime{I <: Integer, T}
    mjd::I
    t0::T
    dt::T
end

mjd(ts::IntegrationTime) = ts.mjd
interval(ts::IntegrationTime) = (ts.t0 - ts.dt / 2) .. (ts.t0 + ts.dt / 2)
Base.in(t::Number, ts::IntegrationTime) = (ts.t0 - ts.dt / 2) ≤ t < (ts.t0 + ts.dt / 2)
Base.isless(t::IntegrationTime, ts::IntegrationTime) = t.t0 < ts.t0
Base.isless(s::Number, t::IntegrationTime) = s < (t.t0 - t.dt / 2)
Base.isless(t::IntegrationTime, s::Number) = (t.t0 + t.dt / 2) < s
Base.Broadcast.broadcastable(ts::IntegrationTime) = Ref(ts)

_center(ts::IntegrationTime) = ts.t0
_region(ts::IntegrationTime) = ts.dt

struct FrequencyChannel{T, I <: Integer}
    central::T
    bandwidth::T
    channel::I
end
channel(fs::FrequencyChannel) = fs.channel
interval(fs::FrequencyChannel) = (fs.central - fs.bandwidth / 2) .. (fs.central + fs.bandwidth / 2)
Base.in(f::Number, fs::FrequencyChannel) = (fs.central - fs.bandwidth / 2) ≤ f < (fs.central + fs.bandwidth / 2)
Base.isless(t::FrequencyChannel, ts::FrequencyChannel) = _center(t) < _center(ts)
Base.isless(s::Number, t::FrequencyChannel) = s < (_center(t) - _region(t) / 2)
Base.isless(t::FrequencyChannel, s::Number) = (_center(t) + _region(t) / 2) < s
Base.Broadcast.broadcastable(fs::FrequencyChannel) = Ref(fs)


_center(fs::FrequencyChannel) = fs.central
_region(fs::FrequencyChannel) = fs.bandwidth


include("site_array.jl")
include("feedrotations.jl")
include("jonesmatrices.jl")
include("priors/priors.jl")
include("instrument_transforms.jl")
include("model.jl")
include("utility.jl")
include("caltable.jl")
