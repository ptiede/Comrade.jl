#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
using Pkg, Distributions
using ChainRulesTestUtils
using Pyehtim
using Optimization
using Comrade
using Test


# Now we need to grab all the subpackages for testing
const GROUP = get(ENV, "GROUP", "ALL")
#using PyCall


# add subpkg
function dev_subpkg(subpkg)
    subpkg_path = joinpath(dirname(@__DIR__), "lib", subpkg)
    Pkg.develop(PackageSpec(path=subpkg_path))
end

# activate subpkg and instantiate it
function activate_subpkg_env(subpkg)
    subpkg_path = joinpath(dirname(@__DIR__), "lib", subpkg)
    Pkg.activate(subpkg_path)
    Pkg.develop(PackageSpec(path=subpkg_path))
    Pkg.instantiate()
end

include(joinpath(@__DIR__, "test_util.jl"))

# Now split depending on what kind of test we are doing
if GROUP == "ALL" || GROUP == "Core"
    Pkg.develop(PackageSpec(url="https://github.com/ptiede/ComradeBase.jl"))
    @testset "Comrade.jl" begin
        include(joinpath(@__DIR__, "Core/core.jl"))
        # include(joinpath(@__DIR__, "ext/comradeahmc.jl"))
        # include(joinpath(@__DIR__, "ext/comradeoptimization.jl"))
        # include(joinpath(@__DIR__, "ext/comradepigeons.jl"))
        # include(joinpath(@__DIR__, "ext/comradedynesty.jl"))
        # include(joinpath(@__DIR__, "ext/comradenested.jl"))
    end
end
