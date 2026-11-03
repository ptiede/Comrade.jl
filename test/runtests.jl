#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
using Pkg, Distributions
using VLBIFiles
using ChainRulesTestUtils
using Optimization
using Comrade
using Test
using LinearAlgebra
using CairoMakie
using Serialization

#using PyCall

include(joinpath(@__DIR__, "test_util.jl"))

# Pkg.develop(PackageSpec(url = "https://github.com/ptiede/ComradeBase.jl"))
@testset "Comrade.jl" begin
    @testset "exported names are defined" begin
        # PolarizedTypes exports some names it does not define.
        upstream = filter(n -> !isdefined(Comrade.PolarizedTypes, n), names(Comrade.PolarizedTypes))
        @test setdiff(filter(n -> !isdefined(Comrade, n), names(Comrade)), upstream) == Symbol[]
    end
    @static if VERSION >= v"1.11"
        @testset "public names are defined and documented" begin
            pub = filter(n -> !Base.isexported(Comrade, n), names(Comrade))
            @test !isempty(pub)
            @test filter(n -> !(isdefined(Comrade, n) && Base.Docs.hasdoc(Comrade, n)), pub) == Symbol[]
        end
    end
    include(joinpath(@__DIR__, "Core/core.jl"))
    include(joinpath(@__DIR__, "ext/comradeahmc.jl"))
    include(joinpath(@__DIR__, "ext/comrademcmcdiagnostictools.jl"))
    include(joinpath(@__DIR__, "ext/comradeoptimization.jl"))
    include(joinpath(@__DIR__, "ext/comradepigeons.jl"))
    include(joinpath(@__DIR__, "ext/comradedynesty.jl"))
    include(joinpath(@__DIR__, "ext/comradenested.jl"))
    include(joinpath(@__DIR__, "ext/comradereactant.jl"))
end
