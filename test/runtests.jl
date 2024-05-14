#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
using SafeTestsets, Pkg, Distributions
using ChainRulesTestUtils
using Pyehtim

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

# Now split depending on what kind of test we are doing
if GROUP == "ALL" || GROUP == "Core"
    dev_subpkg("ComradeOptimization")
    Pkg.develop(PackageSpec(url="https://github.com/ptiede/ComradeBase.jl"))
    @safetestset "CORE Comrade.jl" begin
        using Optimization
        using Comrade
        include(joinpath(@__DIR__, "Core/core.jl"))
        include(joinpath(@__DIR__, "ext/comradeoptimization.jl"))
    end
end
