#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_pycall.jl")
#end
using SafeTestsets, Pkg, Distributions

# Now we need to grab all the subpackages for testing
const GROUP = get(ENV, "GROUP", "ALL")
#using PyCall


# add develop to subpkg
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
if GROUP == "All" || GROUP == "CORE"
    dev_subpkg("ComradeGalactic")
    dev_subpkg("ComradeAHMC")
    @safetestset "CORE Comrade.jl" begin
        include("Core/core.jl")
    end
else
    dev_subpkg(GROUP)
    subpkg_path = joinpath(dirname(@__DIR__), "lib", GROUP)
    Pkg.test(PackageSpec(name=GROUP, path=subpkg_path))
end
