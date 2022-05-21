module ComradeGalactic

using Comrade
using GalacticOptim
using Reexport

@reexport GalacticOptim



function GalacticOptim.OptimizationFunction(post::Comrade.Posterior, args...; kwargs...)
    throw("Transform the posterior first using `asflat` or `ascube`")
end

function GalacticOptim.OptimizationFunction(post::Comrade.TransformedPosterior, args...; kwargs...)
    ℓ(x,p) = -logdensityof(post, x)
    return OptimizationFunction(ℓ, args...; kwargs...)
end

# """
#     solve(prob, opt, args...; transform=nothing, kwargs...)
# Find the optimum of the problem.

# Optional keyword argument `transform``
# will transform the solution to parameter space, and the return will be
# a tuple with the galactic optim solution in the first element and the optimum
# location in model space in the second argument
# """
# function GalacticOptim.solve(prob::OptimizationProblem, opt, transform::Union{Nothing, HypercubeTransform.AbstractHypercubeTransform, HypercubeTransform.TransformVariables.AbstractTransform, HypercubeTransform.NamedFlatTransform}, args...; kwargs...)
#     sol = solve(prob, opt, args...; kwargs...)
#     if isnothing(transform)
#         return sol
#     else
#         p = HypercubeTransform.transform(transform, sol.u)
#         return sol, p
#     end
# end

end
