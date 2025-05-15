module ComradeOptimizationExt

using Comrade

using Optimization
using Distributions
using LinearAlgebra
using HypercubeTransform
using LogDensityProblems
using Optimization: SciMLBase

function Optimization.OptimizationFunction(post::Comrade.TransformedVLBIPosterior, args...; kwargs...)
    ℓ(x, p = post) = -logdensityof(p, x)
    if isnothing(Comrade.admode(post))
        return SciMLBase.OptimizationFunction(ℓ, args...; kwargs...)
    else
        function grad(G, x, p)
            (_, dx) = Comrade.LogDensityProblems.logdensity_and_gradient(post, x)
            dx .*= -1
            G .= dx
            return G
        end
        return SciMLBase.OptimizationFunction(ℓ, args...; grad = grad, kwargs...)
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
 -  `lb`: The lower bounds for the parameters. This is expected to be a NamedTuple of a similar structure as samples from the posterior.
 - `ub` : The upper bounds for the parameters. This is expected to be a NamedTuple of a similar structure as samples from the posterior.
 - `transform`: The transformation of the posterior to use. This can be `ascube` or `asflat` or `nothing`. 
                If `cube` is true then the posterior is transformed to a unit cube. If `asflat` then we transform to support (-∞, ∞).
                If `nothing` then we use heuristic to determine the transformation based on the dimension of the posterior.
 - `kwargs` : Additional keyword arguments passed `Optimization.jl` `solve` method.

"""
function Comrade.comrade_opt(post::VLBIPosterior, opt, args...; initial_params = nothing, lb = nothing, ub = nothing, transform = nothing, kwargs...)

    (!isnothing(lb) & isnothing(ub)) && throw(ArgumentError("You specified lower bounds but not upper bounds. Please specify both or neither."))
    (!isnothing(ub) & isnothing(lb)) && throw(ArgumentError("You specified upper bounds but not lower bounds. Please specify both or neither."))

    tpost = _transform_heuristic(post, opt, lb, transform)
    f = OptimizationFunction(tpost)

    if isnothing(initial_params)
        initial_params = prior_sample(tpost)
    else
        initial_params = Comrade.inverse(tpost, initial_params)
    end


    lbt = nothing
    ubt = nothing

    if tpost.transform isa HypercubeTransform.AbstractHypercubeTransform
        lbt = fill(0.0001, dimension(tpost))
        ubt = fill(0.9999, dimension(tpost))
    end

    if SciMLBase.allowsbounds(opt) && !isnothing(lb)
        lbt = Comrade.inverse(tpost, lb)
        ubt = Comrade.inverse(tpost, ub)
    end


    prob = OptimizationProblem(f, initial_params, tpost; lb = lbt, ub = ubt)
    sol = solve(prob, opt, args...; kwargs...)
    return Comrade.transform(tpost, sol), sol
end

@noinline function _transform_heuristic(post, opt, lb, transform)
    if isnothing(transform)
        # If bounds if no gradient then the cubic transform is the most appropriate
        if SciMLBase.requiresbounds(opt) && !SciMLBase.requiresgradient(opt)
            return ascube(post)
            # If bounds and gradients but the size of the problem is small then we still use the cubic transform
        elseif SciMLBase.requiresbounds(opt) && Comrade.dimension(post) < 100
            return ascube(post)
        elseif SciMLBase.requiresbounds(opt) && !isnothing(lb)
            return asflat(post)
        elseif SciMLBase.requiresbounds(opt) && (isnothing(lb))
            throw(
                ArgumentError(
                    "You are using the $(show(opt)) optimizer which requires bounds, is high dimensional, and wants gradients.\n" *
                        "`asflat` is the best choice here but it requires you to manually specify the bounds.\n" *
                        "Please specify the bounds for the optimization problem."
                )
            )
        elseif !SciMLBase.requiresbounds(opt)
            return asflat(post)
        elseif !SciMLBase.allowsbounds(opt) && (!isnothing(lb))
            @warn "You using the $(show(opt)) optimizer which does not allow bounds but you specified them. We are ignoring them."
            return asflat(post)
        else
            # Default to asflat since we require gradients and now bounds so things is the best option.
            return asflat(post)
        end
    else
        return transform(post)
    end
end

end
