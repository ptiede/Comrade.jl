#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
using Pkg
using PyCall
Pkg.build("PyCall")
println(PyCall.libpython)
using Comrade
using FFTW
using Statistics
using Test

#include("test_util.jl")

@testset "Comrade.jl" begin
    include("observation.jl")
    include("distributions.jl")
    include("models.jl")
end
