# THis is because we are on a headless system
# see https://github.com/jheinen/GR.jl/issues/510
# ENV["GKS_WSTYPE"]="nul"

using Documenter, Pkg
using DocumenterVitepress
using Stoked, StokedBase, AdvancedHMC, Dynesty, Optimization,
      PolarizedTypes
using Pyehtim, VLBISkyModels, InteractiveUtils
using AbstractMCMC, Random, HypercubeTransform
using CairoMakie


Documenter.DocMeta.setdocmeta!(Stoked, :DocTestSetup, :(using Stoked); recursive=true)


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
        repo="https://github.com/ptiede/Stoked.jl",
        devbranch = "main",
        devurl = "dev",
    )


makedocs(;
    modules=[StokedBase, Stoked, Base.get_extension(Stoked, :StokedMakieExt)],
    repo="https://github.com/ptiede/Stoked.jl/blob/{commit}{path}#{line}",
    sitename="Stoked.jl",
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
    repo="github.com/ptiede/Stoked.jl",
    push_preview=true,
    devbranch = "main",
    target = "build"
)
