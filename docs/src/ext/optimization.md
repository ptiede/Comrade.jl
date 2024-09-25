# Optimization Extension

To optimize our posterior, we use the [`Optimization.jl`](https://github.com/SciML/Optimization.jl) package. Optimization provides a global interface to several Julia optimizers. The base call most people should 
look at is [`comrade_opt`](@ref) which serves as the general purpose
optimization algorithm.

To see what optimizers are available and what options are available, please see the `Optimizations.jl` [docs](http://optimization.sciml.ai/dev/).


## Example

```julia
using Comrade
using Optimization
using OptimizationOptimJL
using Enzyme

# Some stuff to create a posterior object
post # of type Comrade.Posterior

xopt, sol = comrade_opt(post, LBFGS(); adtype=Val(:Enzyme))
```