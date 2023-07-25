# THis is because we are on a headless system
# see https://github.com/jheinen/GR.jl/issues/510
# ENV["GKS_WSTYPE"]="nul"


using Documenter
using Pyehtim
using Zygote
using Comrade
using ComradeBase

using Literate
using Pkg

function dev_subpkg(subpkg)
    subpkg_path = joinpath(dirname(@__DIR__), "lib", subpkg)
    Pkg.develop(PackageSpec(path=subpkg_path))
end

# Make sure we are using main branch versions of the packages for the docs
Pkg.develop(PackageSpec(url="https://github.com/ptiede/ComradeBase.jl"))
dev_subpkg("ComradeAHMC")
dev_subpkg("ComradeOptimization")
dev_subpkg("ComradeNested")
dev_subpkg("ComradeDynesty")
dev_subpkg("ComradeAdaptMCMC")

using ComradeAHMC
using ComradeOptimization
using ComradeNested
using ComradeDynesty
using ComradeAdaptMCMC
using PolarizedTypes
using OptimizationBBO
using Glob
using Plots


# Make the examples using Literate
GENERATED = joinpath(@__DIR__, "../", "examples")
OUTDIR = joinpath(@__DIR__, "src", "examples")

SOURCE_FILES = Glob.glob("*.jl", GENERATED)
foreach(fn -> Literate.markdown(fn, OUTDIR, documenter=true), SOURCE_FILES)

MD_FILES = [joinpath("examples", "data.md"),
            joinpath("examples", "geometric_modeling.md"),
            joinpath("examples", "imaging_closures.md"),
            joinpath("examples", "imaging_vis.md"),
            joinpath("examples", "imaging_pol.md"),
            joinpath("examples", "hybrid_imaging.md")
           ]
# joinpath.("examples", replace.(basename.(SOURCE_FILES), ".jl"=>".md"))


makedocs(;
    modules=[ComradeBase, Comrade,
             ComradeOptimization, ComradeAHMC,
             ComradeNested, ComradeDynesty,
             ComradeAdaptMCMC, PolarizedTypes],
    repo="https://github.com/ptiede/Comrade.jl/blob/{commit}{path}#{line}",
    sitename="Comrade.jl",
    pages=Any[
        "Home" => "index.md",
        "benchmarks.md",
        "vlbi_imaging_problem.md",
        "conventions.md",
        "Tutorials" => MD_FILES,
        "Libraries" => [
                        "libs/optimization.md",
                        "libs/ahmc.md",
                        "libs/nested.md",
                        "libs/dynesty.md",
                        "libs/adaptmcmc.md"
                       ],
        "base_api.md",
        "api.md"
    ],
    format = Documenter.HTML(), draft=false
)

deploydocs(;
    repo="github.com/ptiede/Comrade.jl",
    push_preview=false,
    devbranch = "main",
)
