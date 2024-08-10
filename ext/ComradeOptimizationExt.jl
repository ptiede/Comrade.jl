module ComradeOptimizationExt

using Comrade

using Optimization
using Distributions
using LinearAlgebra
using HypercubeTransform


function Optimization.OptimizationFunction(post::Comrade.TransformedVLBIPosterior, args...; kwargs...)
    ℓ(x,p) = -logdensityof(post, x)
    return SciMLBase.OptimizationFunction(ℓ, args...; kwargs...)
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
    comrade_opt(post::VLBIPosterior, opt, adtype=nothing, args...; initial_params=nothing, kwargs...)

Optimize the posterior `post` using the `opt` optimizer.

## Arguments

 - `post` : The posterior to optimize.
 - `opt` : The optimizer to use. This can be any optimizer from `Optimization.jl`.
 - `adtype` : The automatic differentiation type to use. The default is `nothing` which means
    no automatic differentiation is used. To specify to use automatic differentiation
    set `adtype`. For example if you wish to use `Enzyme` set `adtype=Optimization.AutoEnzyme(Enzyme.Reverse)`.
 - `args` : Additional arguments passed to the `Optimization`, `solve` method

## Keyword Arguments

 - `initial_params` : The initial parameters to start the optimization. If `nothing` then
    the initial parameters are sampled from the prior. If not `nothing` then the initial
    parameters are transformed to the transformed space.
 - `kwargs` : Additional keyword arguments passed `Optimization.jl` `solve` method.

"""
function Comrade.comrade_opt(post::VLBIPosterior, opt, adtype=nothing, args...; initial_params=nothing, kwargs...)
    if isnothing(adtype)
        adtype = Optimization.SciMLBase.NoAD()
        tpost = ascube(post)
    else
        tpost = asflat(post)
    end

    f = OptimizationFunction(tpost, adtype)

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

    prob = OptimizationProblem(f, initial_params, nothing; lb, ub)
    sol = solve(prob, opt, args...; kwargs...)
    return transform(tpost, sol), sol
end

end
