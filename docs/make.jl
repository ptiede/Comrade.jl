using Documenter
using Comrade
using ComradeBase

using Literate

using Pkg
Pkg.develop(path=joinpath(dirname(@__DIR__), "."))
Pkg.develop(path=joinpath(dirname(@__DIR__), "lib", "ComradeAHMC"))
Pkg.develop(path=joinpath(dirname(@__DIR__), "lib", "ComradeOptimization"))
Pkg.develop(path=joinpath(dirname(@__DIR__), "lib", "ComradeAdaptMCMC"))
Pkg.develop(path=joinpath(dirname(@__DIR__), "lib", "ComradeDynesty"))
Pkg.develop(path=joinpath(dirname(@__DIR__), "lib", "ComradeNested"))

using ComradeAHMC
using ComradeOptimization
using ComradeNested
using ComradeDynesty
using ComradeAdaptMCMC
using OptimizationBBO
using Glob
using Plots


# Make the examples using Literate
GENERATED = joinpath(@__DIR__, "src", "examples")
SOURCE_FILES = Glob.glob("*.jl", GENERATED)
foreach(fn -> Literate.markdown(fn, GENERATED), SOURCE_FILES)



format = Documenter.HTML(edit_link = "source",
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = String[])


ENV["GKSwstype"] = "nul" # needed for the GR backend on headless servers
gr()



makedocs(;
    modules=[ComradeBase, Comrade,
             ComradeOptimization, ComradeAHMC,
             ComradeNested, ComradeDynesty,
             ComradeAdaptMCMC],
    authors="Paul Tiede",
    repo="https://github.com/ptiede/Comrade.jl/blob/{commit}{path}#L{line}",
    sitename="Comrade.jl",
    format=format,
    pages=Any[
        "Home" => "index.md",
        "benchmarks.md",
        "vlbi_imaging_problem.md",
        "Tutorials" => [
                       "examples/data.md",
                       "examples/black_hole_image.md",
                       "examples/nonanalytic.md"
                       ],
        "Libraries" => [
                        "libs/optimization.md",
                        "libs/ahmc.md",
                        "libs/nested.md",
                        "libs/dynesty.md",
                        "libs/adaptmcmc.md"
                       ],
        "interface.md",
        "base_api.md",
        "api.md"
    ],
)

deploydocs(;
    repo="github.com/ptiede/Comrade.jl",
    push_preview=true,
    forcepush=true,
    devbranch = "main",
)
