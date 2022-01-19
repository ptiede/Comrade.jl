#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
using Pkg
#ENV["PYTHON"]="/opt/hostedtoolcache/Python/3.10.1/x64/lib/libpython3.10.so.1.0"
#using PyCall
#run(`$(PyCall.python) -m pip install --upgrade pip setuptools`)
#run(`$(PyCall.python) -m pip install ehtim`)
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
