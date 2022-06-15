module ComradeGalactic

using Comrade
using Reexport
using Distributions
using ForwardDiff
using LinearAlgebra

@reexport using GalacticOptim

export laplace


function GalacticOptim.OptimizationFunction(post::Comrade.Posterior, args...; kwargs...)
    throw("Transform the posterior first using `asflat` or `ascube`")
end

function GalacticOptim.OptimizationFunction(post::Comrade.TransformedPosterior, args...; kwargs...)
    ℓ(x,p) = -logdensityof(post, x)
    return OptimizationFunction(ℓ, args...; kwargs...)
end

"""
    `laplace(prob, opt, args...; kwargs...)`
Compute the Laplace or Quadratic approximation to the prob or posterior.
The `args` and `kwargs` are passed the the GalacticOptim.solve function.

Note the quadratic approximation is in the space of the transformed posterior
not the usual parameter space. This is better for constrained problems where
we may run up against a boundary.
"""
function laplace(prob::OptimizationProblem, opt, args...; kwargs...)
    sol = solve(prob, opt, args...; kwargs...)
    f = Base.Fix2(prob.f, nothing)
    J = ForwardDiff.hessian(f, sol)
    h = J*sol
    return MvNormalCanon(h, Symmetric(J))
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
