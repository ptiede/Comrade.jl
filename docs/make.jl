# THis is because we are on a headless system
# see https://github.com/jheinen/GR.jl/issues/510
# ENV["GKS_WSTYPE"]="nul"

using Documenter, Pkg
using DocumenterVitepress
using Comrade, ComradeBase, AdvancedHMC, Dynesty, NestedSamplers, Optimization,
      PolarizedTypes
using Pyehtim, VLBISkyModels, InteractiveUtils
using AbstractMCMC, Random, HypercubeTransform


deployconfig = Documenter.auto_detect_deploy_system()
Documenter.post_status(deployconfig; type="pending", repo="github.com/ptiede/Comrade.jl.git")

TUTORIALS = [
        "tutorials/ClosureImaging.md",
        "tutorials/GeometricModeling.md",
        "tutorials/HybridImaging.md",
        "tutorials/LoadingData.md",
        "tutorials/PolarizedImaging.md",
        "tutorials/StokesIImaging.md",
     ]

makedocs(;
    modules=[ComradeBase, Comrade],
    repo="https://github.com/ptiede/Comrade.jl/blob/{commit}{path}#{line}",
    sitename="Comrade.jl",
    format = MarkdownVitepress(
        repo="https://github.com/ptiede/Comrade.jl",
        devurl = "dev",
        devbranch = "main",
    ),
    pages=Any[
        "Home" => "index.md",
        "benchmarks.md",
        "vlbi_imaging_problem.md",
        "conventions.md",
        "Tutorials" => TUTORIALS,
        "Extensions" => [
                        "ext/optimization.md",
                        "ext/ahmc.md",
                        "ext/nested.md",
                        "ext/dynesty.md",
                       ],
        "base_api.md",
        "api.md"
    ],
    draft = false,
    source = "src",
    build  = "build"
)

deploydocs(;
    repo="github.com/ptiede/Comrade.jl",
    push_preview=true,
    devbranch = "main",
    target = "build"
)
