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
GENERATED = joinpath(@__DIR__, "../", "examples")
OUTDIR = joinpath(@__DIR__, "src", "examples")

SOURCE_FILES = Glob.glob("*.jl", GENERATED)
println(SOURCE_FILES)
foreach(fn -> Literate.markdown(fn, OUTDIR, documenter=true), SOURCE_FILES)

MD_FILES = joinpath.("examples", replace.(basename.(SOURCE_FILES), ".jl"=>".md"))

deployconfig = Documenter.auto_detect_deploy_system()
Documenter.post_status(deployconfig; type="pending", repo="github.com/avik-pal/Lux.jl.git")



makedocs(;
    modules=[ComradeBase, Comrade,
             ComradeOptimization, ComradeAHMC,
             ComradeNested, ComradeDynesty,
             ComradeAdaptMCMC],
    sitename="Comrade.jl",
    pages=Any[
        "Home" => "index.md",
        "benchmarks.md",
        "vlbi_imaging_problem.md",
        "Tutorials" => MD_FILES,
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
    build = joinpath(@__DIR__, "docs"), draft = true,
    format = Documenter.HTML()
)

deploydocs(;
    repo="github.com/ptiede/Comrade.jl",
    push_preview=true,
    devbranch = "main",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", devurl => devurl],
)
