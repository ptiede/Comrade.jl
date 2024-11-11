using MacroTools
using NamedTupleTools


# Idea
#  1. capture all tilde statements and convert them to a NamedTupleDist
#  2. capture all equals signs and convert make a function with a gensym name
#  3. capture all arguments and make a named tuple with them as metadata


# What I want to do
# @sky function gauss(g; a, b, c)
#      ax ~ Uniform (0.0, a)
#      bx ~ Uniform (0.0, b)
#      cx ~ Uniform (0.0, c)
#      m = modify(Gaussian(), Stretch(ax, bx), Renormlize(cx))
#      return m
#end 
# 
# Parses into
# function __gaussskymod__(θ, metadata)
#        (; ax, bx, cx) = θ
#        (; a, b, c, g) = metadata
#        m = modify(Gaussian(), Stretch(ax, bx), Renormlize(cx))
#        return m
# end
# gauss = function(g; a, b, c)
#     metadata = (; a, b, c, g=g)
#     prior = (
#            ax = Uniform(0.0, a),
#            bx = Uniform(0.0, b),
#            cx = Uniform(0.0, c)
#           )
#    
#    return SkyModel(__gaussskymod__, prior, g; metadata=metadata)
#end

macro sky(expr)
    # println(expr)
    def = splitdef(expr)
    name = def[:name]
    args = def[:args]
    kwargs = def[:kwargs]
    body = def[:body]

    @info kwargs
    kwvar = map(x->x.args[1], kwargs)
    kwval = map(x->x.args[2], kwargs)

    kwdef = NamedTuple{Tuple(kwvar)}(Tuple(kwval))
    gexpr = Expr(:kw, :g, esc(g))
    metadata = build_metadata([kwargs...,gexpr])
    prior    = build_prior(metadata, body)
    # skym     = build_sky(args, kwargs, body)

    model_def = MacroTools.combinedef(skymodel)

    return quote
        $name = function($g; $kwdef...)
            metadata = (; kwargs..., g=$g)
        end
    end
end

function build_metadata(kwargs)
    isempty(kwargs) && return nothing
    return Expr(:tuple, Expr(:parameters, kwargs...,))
end

function build_prior(metadata, body)
    names = []
    dists = []
    postwalk(x->@capture(x, T_ ~ D_) ? (push!(names, T);push!(dists, D)) : x, ex)
    return Expr(:call, :(HypercubeTransform.NamedDist), Expr(:parameters, ))
    return :(NamedDist(NamedTuple{Tuple($names)}($dists)))
end
