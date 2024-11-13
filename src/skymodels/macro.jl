using MacroTools
using NamedTupleTools
using HypercubeTransform


# Idea
#  1. capture all tilde statements and convert them to a NamedTupleDist
#  2. capture all equals signs and convert make a function with a gensym name
#  3. capture all arguments and make a named tuple with them as metadata


# What I want to do
# @sky function gauss(g; a=1.0, b=1.0, c=2.0)
#      ax ~ Uniform (0.0, a)
#      bx ~ Uniform (0.0, b)
#      cx ~ Uniform (0.0, c)
#      m = modify(Gaussian(), Stretch(ax, bx), Renormlize(cx))
#      return m
#end 
# 
# Parses into
# function __gaussskymod__(θ, (; a = 1.0, b = 1.0, c = 2.0, grid=g))
#        (; ax, bx, cx) = θ
#        (; a, b, c, g) = metadata
#        m = modify(Gaussian(), Stretch(ax, bx), Renormlize(cx))
#        return m
# end
# 
# function __gaussskyprior__(;a = 1.0, b = 1.0, c = 2.0, grid=g)
#     ax = Uniform(0.0, a)
#     bx = Uniform(0.0, b)
#     cx = Uniform(0.0, c)
#     return HypercubeTransform.NamedDist((;ax, bx, cx))
#end


# gauss = function(g; a, b, c)
#     metadata = (; a, b, c, grid=g)
#    return SkyModel(__gaussskymod__, __gaussskyprior__(;metadata...), g; metadata=metadata)
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
    prior     = build_prior(name, metadata, body)
    skymodel  = build_sky(name, args, kwargs, body)
    @info metadata
    @info prior
    return prior

    skyprf = Symbol("__", name, "skyprior__")

    model_def = MacroTools.combinedef(skymodel)

    return quote
        $skyprf
        $skym
    end
end

function build_metadata(kwargs)
    isempty(kwargs) && return nothing
    return Expr(:tuple, Expr(:parameters, kwargs...,))
end

function build_prior(metadata, body)
    names = []
    dists = []
    MacroTools.postwalk(x->@capture(x, T_ ~ D_) ? (push!(names, T);push!(dists, D)) : x, body)
    body = 
    return Expr(:call, :(HypercubeTransform.NamedDist), Expr(:parameters, names, dists))
    # return :(HypercubeTransform.NamedDist(NamedTuple{Tuple($names)}($dists)))
end
