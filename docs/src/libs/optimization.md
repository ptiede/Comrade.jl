# ComradeOptimization

To optimize our posterior we use the [`Optimization.jl`](https://github.com/SciML/Optimization.jl) package. This provides a global interface to several Julia optimizers. The `Comrade` wrapper for `Optimization.jl` is very thin. The only difference addition is that `Comrade` has provided a method:

```julia
OptimizationFunction(::TransformedPosterior, args...; kwargs...)
```

meaning we can pass it a posterior object and it will set up the `OptimizationFunction` for us. **Note** that we only specify this for 
a transformed version of the posterior. This is because `Optimization.jl` requires a flattened version of the posterior. Additionally different optimizers may prefer different parameter transformations. For examples if we use `OptimizationBBO` using [`ascube`](@ref) is a good choice since it needs a compact region to search over, and `ascube` convert our parameter space to the unit hypercube. On the other hand gradient based optimizers work best without bounds, so a better choice there would be the [`asflat`](@ref) transformation.

To see what optimizers are available and what options are available please see the `Optimizations.jl` [docs](http://optimization.sciml.ai/dev/).

## Example

```julia
using Comrade
using ComradeOptimization
using OptimizationOptimJL

# Some stuff to create a posterior object
post # of type Comrade.Posterior

# Create a optimization problem using ForwardDiff as the backend
fflat = OptimizationProblem(asflat(post), Optimization.AutoForwardDiff())

# create the problem from a random point in the prior, nothing is b/c there are no additional arugments to our function.
prob = OptimizationProblem(fflat, prior_sample(asflat(post)), nothing)

# Now solve! Here we use LBFGS
sol = solve(prob, LBFGS(); g_tol=1e-2)
```


## API

```@meta
CurrentModule = ComradeOptimization
```

```@autodocs
Modules = [ComradeOptimization]
Order   = [:function, :type]
```
