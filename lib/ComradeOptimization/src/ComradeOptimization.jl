module ComradeOptimization

using Comrade
using Reexport
using Distributions
using ForwardDiff
using LinearAlgebra
import SciMLBase

@reexport using Optimization

export laplace


function SciMLBase.OptimizationFunction(post::Comrade.Posterior, args...; kwargs...)
    throw("Transform the posterior first using `asflat` or `ascube`")
end

"""
    SciMLBase.OptimizationFunction(post::Posterior, args...; kwargs...)

Constructs a `OptimizationFunction` from a `Comrade.TransformedPosterior` object.
Note that a user must **transform the posterior first**. This is so we know which
space is most amenable to optimization.
"""
function SciMLBase.OptimizationFunction(post::Comrade.TransformedPosterior, args...; kwargs...)
    ℓ(x,p) = -logdensityof(post, x)
    return SciMLBase.OptimizationFunction(ℓ, args...; kwargs...)
end

"""
    laplace(prob, opt, args...; kwargs...)

Compute the Laplace or Quadratic approximation to the prob or posterior.
The `args` and `kwargs` are passed the the SciMLBase.solve function.
This will return a `Distributions.MvNormal` object that approximates
the posterior in the transformed space.

Note the quadratic approximation is in the space of the transformed posterior
not the usual parameter space. This is better for constrained problems where
we may run up against a boundary.
"""
function laplace(prob::SciMLBase.OptimizationProblem, opt, args...; kwargs...)
    sol = solve(prob, opt, args...; kwargs...)
    f = Base.Fix2(prob.f, nothing)
    J = ForwardDiff.hessian(f, sol)
    @. J += 1e-5 # add offset to help with positive definitness
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
