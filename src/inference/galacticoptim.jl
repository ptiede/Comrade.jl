using .GalacticOptim

function GalacticOptim.OptimizationFunction(post::Posterior, args...; kwargs...)
    tpost = asflat(post)
    ℓ(x,p) = -logdensity(tpost, x)
    return OptimizationFunction(ℓ, args...; kwargs...), tpost.transform
end

function GalacticOptim.OptimizationFunction(post::TransformedPosterior, args...; kwargs...)
    ℓ(x,p) = -logdensity(post, x)
    return OptimizationFunction(ℓ, args...; kwargs...), post.transform
end

"""
    solve(prob, opt, args...; transform=nothing, kwargs...)
Find the optimum of the problem.

Optional keyword argument `transform``
will transform the solution to parameter space, and the return will be
a tuple with the galactic optim solution in the first element and the optimum
location in model space in the second argument
"""
function GalacticOptim.solve(prob::OptimizationProblem, opt, transform, args...; kwargs...)
    sol = solve(prob, opt, args...; kwargs...)
    if isnothing(transform)
        return sol
    else
        p = HypercubeTransform.transform(transform, sol.u)
        return sol, p
    end
end
