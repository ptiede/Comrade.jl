# THis is because we are on a headless system
# see https://github.com/jheinen/GR.jl/issues/510
# ENV["GKS_WSTYPE"]="nul"

using Documenter, Pkg
using Comrade, ComradeBase,
     ComradeAHMC,
     ComradeOptimization,
     ComradeNested,
     ComradeDynesty,
     ComradeAdaptMCMC

using Pyehtim, VLBISkyModels, InteractiveUtils

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
    modules=[ComradeBase, Comrade,
             ComradeOptimization, ComradeAHMC,
             ComradeNested, ComradeDynesty,
             ComradeAdaptMCMC],
    # repo="https://github.com/ptiede/Comrade.jl/blob/{commit}{path}#{line}",
    sitename="Comrade.jl",
    format = Documenter.HTML(;size_threshold_ignore=TUTORIALS
    ),
    pages=Any[
        "Home" => "index.md",
        "benchmarks.md",
        "vlbi_imaging_problem.md",
        "conventions.md",
        #"Tutorials" => TUTORIALS,
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
)

deploydocs(;
    repo="github.com/ptiede/Comrade.jl",
    push_preview=false,
    devbranch = "main",
)
