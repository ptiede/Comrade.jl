using ROSE
using Documenter

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
)
