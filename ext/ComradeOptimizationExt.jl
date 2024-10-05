module ComradeOptimizationExt

using Comrade

using Optimization
using Distributions
using LinearAlgebra
using HypercubeTransform
using LogDensityProblems

function Optimization.OptimizationFunction(post::Comrade.TransformedVLBIPosterior, args...; kwargs...)
    ℓ(x,p=post) = -logdensityof(p, x)
    if isnothing(Comrade.admode(post))
        return SciMLBase.OptimizationFunction(ℓ, args...; kwargs...)
    else
        function grad(G, x, p)
            (_, dx) = LogDensityProblems.logdensity_and_gradient(post, x)
            dx .*= -1
            G .= dx
            return G
        end
        return SciMLBase.OptimizationFunction(ℓ, args...; grad=grad, kwargs...)
    end
end

# """
#     comrade_laplace(prob, opt, args...; kwargs...)

# Compute the Laplace or Quadratic approximation to the prob or posterior.
# The `args` and `kwargs` are passed the the SciMLBase.solve function.
# This will return a `Distributions.MvNormal` object that approximates
# the posterior in the transformed space.

# Note the quadratic approximation is in the space of the transformed posterior
# not the usual parameter space. This is better for constrained problems where
# we may run up against a boundary.
# """
# function Comrade.comrade_laplace(post::VLBIPosterior, opt, adtype=Val(:ForwardDiff))
#     sol = solve(prob, opt, args...; kwargs...)
#     f = Base.Fix2(prob.f, nothing)
#     J = ForwardDiff.hessian(f, sol)
#     @. J += 1e-5 # add offset to help with positive definitness
#     h = J*sol
#     return MvNormalCanon(h, Symmetric(J))
# end


"""
    comrade_opt(post::VLBIPosterior, opt, args...; initial_params=nothing, kwargs...)

Optimize the posterior `post` using the `opt` optimizer.

!!! warning
    To use use a gradient optimizer with AD, `VLBIPosterior` must be created with a specific `admode` specified.
    The `admode` can be a union of `Nothing` and `<:EnzymeCore.Mode` types. We recommend
    using `Enzyme.set_runtime_activity(Enzyme.Reverse)`.


## Arguments

 - `post` : The posterior to optimize.
 - `opt` : The optimizer to use. This can be any optimizer from `Optimization.jl`.
 - `args` : Additional arguments passed to the `Optimization`, `solve` method

## Keyword Arguments

 - `initial_params` : The initial parameters to start the optimization. If `nothing` then
    the initial parameters are sampled from the prior. If not `nothing` then the initial
    parameters are transformed to the transformed space.
 - `kwargs` : Additional keyword arguments passed `Optimization.jl` `solve` method.

"""
function Comrade.comrade_opt(post::VLBIPosterior, opt, args...; initial_params=nothing, kwargs...)
    if isnothing(Comrade.admode(post))
        tpost = ascube(post)
    else
        tpost = asflat(post)
    end

    f = OptimizationFunction(tpost)

    if isnothing(initial_params)
        initial_params = prior_sample(tpost)
    else
        initial_params = Comrade.inverse(tpost, initial_params)
    end

    lb = nothing
    ub = nothing
    if tpost.transform isa HypercubeTransform.AbstractHypercubeTransform
        lb=fill(0.0001, dimension(tpost))
        ub = fill(0.9999, dimension(tpost))
    end

    prob = OptimizationProblem(f, initial_params, tpost; lb, ub)
    sol = solve(prob, opt, args...; kwargs...)
    return transform(tpost, sol), sol
end

end
