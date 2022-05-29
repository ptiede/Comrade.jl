#println(PyCall.libpython)
using Comrade
using FFTW
using Statistics
using Test
using Plots

include(joinpath(@__DIR__, "observation.jl"))
include(joinpath(@__DIR__, "distributions.jl"))
include(joinpath(@__DIR__, "models.jl"))
include(joinpath(@__DIR__, "bayes.jl"))
