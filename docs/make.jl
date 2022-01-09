using Documenter

#if haskey(ENV, "GITHUB_ACTIONS")
#    ENV["JULIA_DEBUG"] = "Documenter"
#end

#Documenter.post_status(; type="pending", repo="github.com/ptiede/ROSE.jl.git")
using ROSE


makedocs(;
    modules=[ROSE],
    authors="Paul Tiede",
    repo="https://github.com/ptiede/ROSE.jl/blob/{commit}{path}#L{line}",
    sitename="ROSE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ptiede.github.io/ROSE.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ptiede/ROSE.jl",
    push_preview=true,
    forcepush=true,
    devbranch = "main",
    devurl = "latest",
)
