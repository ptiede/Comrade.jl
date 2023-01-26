using DelimitedFiles
using AstroTime: modified_julian

export uvpositions, stations, data, arrayconfig,
       scans, bandwidth, telescope_array
       getuv, baselines, scantable


abstract type Observation{T} end

include(joinpath(@__DIR__, "io.jl"))
include(joinpath(@__DIR__, "ehtim.jl"))
