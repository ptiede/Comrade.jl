using MacroTools
using NamedTupleTools


# Idea
#  1. capture all tilde statements and convert them to a NamedTupleDist
#  2. capture all equals signs and convert make a function with a gensym name
#  3. capture all arguments and make a named tuple with them as metadata


macro sky(expr)
    # println(expr)
    def = splitdef(expr)
    name = def[:name]
    args = def[:args]
    kwargs = def[:kwargs]
    body = def[:body]

    g = args[1]
    gexpr = Expr(:kw, g, esc(g))
    metadata = build_metadata([kwargs...,gexpr])
    prior    = build_prior(metadata, body)
    # skym     = build_sky(args, kwargs, body)
    return quote
        $metadata
        # $(esc(name)) = SkyModel($skym, $prior, $g; metadata=$metadata)
        # $name = SkyModel($skym, $prior, $g; metadata=$metadata)
    end
end

function build_metadata(kwargs)
    isempty(kwargs) && return nothing
    return Expr(:tuple, Expr(:parameters, kwargs..., (:)))
end

function build_prior(metadata, body)

end
