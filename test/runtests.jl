#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
using Pkg, Distributions
using ChainRulesTestUtils
using Pyehtim
using Optimization
using Comrade
using Test

#using PyCall

include(joinpath(@__DIR__, "test_util.jl"))

Pkg.develop(PackageSpec(url="https://github.com/ptiede/ComradeBase.jl"))
@testset "Comrade.jl" begin
    include(joinpath(@__DIR__, "Core/core.jl"))
    include(joinpath(@__DIR__, "ext/comradeahmc.jl"))
    include(joinpath(@__DIR__, "ext/comradeoptimization.jl"))
    include(joinpath(@__DIR__, "ext/comradepigeons.jl"))
    include(joinpath(@__DIR__, "ext/comradedynesty.jl"))
    include(joinpath(@__DIR__, "ext/comradenested.jl"))
end
