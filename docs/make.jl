using Documenter

#if haskey(ENV, "GITHUB_ACTIONS")
#    ENV["JULIA_DEBUG"] = "Documenter"
#end

#Documenter.post_status(; type="pending", repo="github.com/ptiede/ROSEx.jl.git")
using ROSEx


makedocs(;
    modules=[ROSEx],
    authors="Paul Tiede",
    repo="https://github.com/ptiede/ROSEx.jl/blob/{commit}{path}#L{line}",
    sitename="ROSEx.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ptiede.github.io/ROSEx.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ptiede/ROSEx.jl",
    push_preview=true,
    forcepush=true
)
