using Documenter
using Comrade

using Literate

using AdvancedHMC
using GalacticOptim
using Glob
using Pathfinder
using Plots

# branch = "master"

# demos, demos_cb, demo_assets = makedemos("src/examples"; branch = branch)

# assets = []
# isnothing(demo_assets) || push!(assets, demo_assets)

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
    modules=[Comrade],
    authors="Paul Tiede",
    repo="https://github.com/ptiede/Comrade.jl/blob/{commit}{path}#{line}",
    sitename="Comrade.jl",
    format=format,
    pages=Any[
        "Home" => "index.md",
        "vlbi_imaging_problem.md",
        "Tutorials" => [
                       "examples/data.md",
                       "examples/black_hole_image.md"
                       ],
        "interface.md",
        "api.md"
    ],
)

deploydocs(;
    repo="github.com/ptiede/Comrade.jl",
    push_preview=true,
    forcepush=true,
    devbranch = "main",
)
