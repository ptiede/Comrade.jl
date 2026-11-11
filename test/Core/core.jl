#println(PyCall.libpython)
using Comrade
using FFTW
using Statistics
using Test
using VLBIImagePriors

include(joinpath(@__DIR__, "plots.jl"))
include(joinpath(@__DIR__, "selfcal.jl"))
include(joinpath(@__DIR__, "gauss_markov.jl"))
include(joinpath(@__DIR__, "models.jl"))
include(joinpath(@__DIR__, "imgnormal.jl"))
include(joinpath(@__DIR__, "bayes.jl"))
include(joinpath(@__DIR__, "fields.jl"))
include(joinpath(@__DIR__, "posterior_data.jl"))
include(joinpath(@__DIR__, "posteriorsamples.jl"))
include(joinpath(@__DIR__, "moves.jl"))
include(joinpath(@__DIR__, "geweke.jl"))
include(joinpath(@__DIR__, "macro.jl"))
