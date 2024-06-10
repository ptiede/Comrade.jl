using Pkg; Pkg.activate(@__DIR__)
using Literate

preprocess(path, str) = replace(str, "__DIR = @__DIR__" => "__DIR = \"$(dirname(path))\"")

get_example_path(p) = joinpath(@__DIR__, "..", "examples", p)
OUTPUT = joinpath(@__DIR__, "src", "tutorials")


TUTORIALS = [
        "beginner/LoadingData/main.jl",
        # "beginner/GeometricModeling/main.jl",
        # "intermediate/ClosureImaging/main.jl",
        # "intermediate/StokesIImaging/main.jl",
        # "intermediate/PolarizedImaging/main.jl",
        # "advanced/HybridImaging/main.jl",
        ]

withenv("JULIA_DEBUG"=>"Literate") do
    for (d, paths) in (("", TUTORIALS),),
        (i,p) in enumerate(paths)
        println(p)
        name = "$((rsplit(p, "/")[2]))"
        d    = "$((rsplit(p, "/")[1]))"
        p_ = get_example_path(p)
        @info p_
        jl_expr = "using Literate;"*
                  "preprocess(path, str) = replace(str, \"__DIR = @__DIR__\" => \"__DIR = \\\"\$(dirname(path))\\\"\");"*
                  "Literate.markdown(\"$(p_)\", \"$(joinpath(OUTPUT, d))\";"*
                  "name=\"$name\", execute=true, flavor=Literate.DocumenterFlavor(),"*
                  "preprocess=Base.Fix1(preprocess, \"$(p_)\"))"
        cm = `julia --project=$(@__DIR__) -e $(jl_expr)`
        run(cm)

    end
end
