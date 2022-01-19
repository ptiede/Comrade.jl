#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
using Pkg
Pkg.build("PyCall")
using PyCall
run(`$(PyCall.python) -m pip install --upgrade pip setuptools`)
run(`$(PyCall.python) -m pip install ehtim`)
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
