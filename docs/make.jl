using Documenter

#if haskey(ENV, "GITHUB_ACTIONS")
#    ENV["JULIA_DEBUG"] = "Documenter"
#end

#Documenter.post_status(; type="pending", repo="github.com/ptiede/Comrade.jl.git")
using Comrade


makedocs(;
    modules=[Comrade],
    authors="Paul Tiede",
    repo="https://github.com/ptiede/Comrade.jl/blob/{commit}{path}#L{line}",
    sitename="Comrade.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ptiede.github.io/Comrade.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "vlbi_imaging_problem.md",
        "interface.md",
        "api.md"
    ],
)

deploydocs(;
    repo="github.com/ptiede/Comrade.jl",
    push_preview=true,
    forcepush=true,
    devbranch = "main",
    devurl = "latest",
)
