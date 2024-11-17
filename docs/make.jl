# THis is because we are on a headless system
# see https://github.com/jheinen/GR.jl/issues/510
# ENV["GKS_WSTYPE"]="nul"

using Documenter, Pkg
using DocumenterVitepress
using Comrade, ComradeBase, AdvancedHMC, Dynesty, Optimization,
      PolarizedTypes
using Pyehtim, VLBISkyModels, InteractiveUtils
using AbstractMCMC, Random, HypercubeTransform


Documenter.DocMeta.setdocmeta!(Krang, :DocTestSetup, :(using Comrade); recursive=true)


TUTORIALS = [
        "Overview" => "tutorials/index.md",
        "Beginner" =>[
            "tutorials/beginner/LoadingData.md",
            "tutorials/beginner/GeometricModeling.md"
            ],
        "Intermediate" => [
           "tutorials/intermediate/ClosureImaging.md",
           "tutorials/intermediate/StokesIImaging.md",
           "tutorials/intermediate/PolarizedImaging.md"
           ],
        "Advanced" => [
           "tutorials/advanced/HybridImaging.md",
           ]
     ]

format = MarkdownVitepress(
        repo="https://github.com/ptiede/Comrade.jl",
        devbranch = "main",
        devurl = "dev",
    )


makedocs(;
    modules=[ComradeBase, Comrade],
    repo="https://github.com/ptiede/Comrade.jl/blob/{commit}{path}#{line}",
    sitename="Comrade.jl",
    format = format,
    draft  = false,
    source = "src",
    build  = "build",
    pages=Any[
        "Home" => "index.md",
        "introduction.md",
        "benchmarks.md",
        "vlbi_imaging_problem.md",
        "conventions.md",
        "Tutorials" => TUTORIALS,
        "Extensions" => [
                        "ext/optimization.md",
                        "ext/ahmc.md",
                        "ext/dynesty.md",
                        "ext/pigeons.md"
                       ],
        "base_api.md",
        "api.md"
    ],
)

deploydocs(;
    repo="github.com/ptiede/Comrade.jl",
    push_preview=true,
    devbranch = "main",
    target = "build"
)
