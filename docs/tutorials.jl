using Pkg
Pkg.activate(@__DIR__)
using Literate

const TUTORIALS = [
    "beginner/LoadingData/main.jl",
    "beginner/GeometricModeling/main.jl",
    "intermediate/ClosureImaging/main.jl",
    "intermediate/StokesIImaging/main.jl",
    "intermediate/PolarizedImaging/main.jl",
    "advanced/Hibi/main.jl",
    "advanced/FitPS/main.jl",
]

const OUTPUT = joinpath(@__DIR__, "src", "tutorials")
get_example_path(p) = joinpath(@__DIR__, "..", "examples", p)

"""
    build_tutorial(rel_path)

Render a single tutorial (e.g. `"beginner/LoadingData/main.jl"`) to markdown
under `docs/src/tutorials/<area>/<name>.md`. Runs in-process so the caller
already has the docs environment loaded.
"""
function build_tutorial(rel_path)
    p = get_example_path(rel_path)
    parts = rsplit(rel_path, "/")
    area = String(parts[1])
    name = String(parts[2])
    preprocess(s) = replace(s, "__DIR = @__DIR__" => "__DIR = \"$(dirname(p))\"")
    return Literate.markdown(
        p, joinpath(OUTPUT, area);
        name = name,
        execute = true,
        flavor = Literate.DocumenterFlavor(),
        preprocess = preprocess,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    withenv("JULIA_DEBUG" => "Literate") do
        for p in TUTORIALS
            @info "Building $p"
            try
                build_tutorial(p)
            catch e
                @warn "Failed to build $p" exception = e
            end
        end
    end
end
