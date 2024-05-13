using LinearAlgebra
using SparseArrays
import Distributions
using Statistics
using PrettyTables

struct IntegrationTime{T}
    mjd::Int
    t0::T
    dt::T
end

Base.in(t::Number, ts::IntegrationTime) = (ts.t0 - ts.dt/2) ≤ t < (ts.t0 + ts.dt)
Base.isless(t::IntegrationTime, ts::IntegrationTime) = t.t0 < ts.t0
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
include("model.jl")
include("utility.jl")
include("caltable.jl")
