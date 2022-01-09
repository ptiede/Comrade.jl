#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
#using Pkg
#Pkg.build("PyCall")
using ROSE
using Statistics
using Test

include("test_util.jl")

@testset "ROSE.jl" begin

    #include("observation.jl")
    include("distributions.jl")
end
