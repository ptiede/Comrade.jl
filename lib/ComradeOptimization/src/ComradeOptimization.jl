module ComradeOptimization

using Comrade: VLBIPosterior, asflat, ascube, transform
using Reexport
using Distributions
using ForwardDiff
using LinearAlgebra
import SciMLBase

@reexport using Optimization


function __init__()
    @warn "ComradeOptimization is deprecated. Optimization.jl is now an extension."
end

# """
#     solve(prob, opt, args...; transform=nothing, kwargs...)
# Find the optimum of the problem.

# Optional keyword argument `transform``
# will transform the solution to parameter space, and the return will be
# a tuple with the galactic optim solution in the first element and the optimum
# location in model space in the second argument
# """
# function SciMLBase.solve(prob::OptimizationProblem, opt, transform::Union{Nothing, HypercubeTransform.AbstractHypercubeTransform, HypercubeTransform.TransformVariables.AbstractTransform, HypercubeTransform.NamedFlatTransform}, args...; kwargs...)
#     sol = solve(prob, opt, args...; kwargs...)
#     if isnothing(transform)
#         return sol
#     else
#         p = HypercubeTransform.transform(transform, sol.u)
#         return sol, p
#     end
# end

end
