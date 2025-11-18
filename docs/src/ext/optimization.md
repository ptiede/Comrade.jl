# Optimization Extension

To optimize our posterior, we use the [`Optimization.jl`](https://github.com/SciML/Optimization.jl) package. Optimization provides a global interface to several Julia optimizers. The base call most people should 
look at is [`comrade_opt`](@ref) which serves as the general purpose
optimization algorithm.

To see what optimizers are available and what options are available, please see the `Optimizations.jl` [docs](http://optimization.sciml.ai/dev/).


!!! warning
    To use use a gradient optimizer with AD, `VLBIPosterior` must be created with a specific `admode` specified.
    The `admode` can be a union of `Nothing` and `<:EnzymeCore.Mode` types. We recommend
    using `Enzyme.set_runtime_activity(Enzyme.Reverse)`


## Example

```julia
using Comrade
using Optimization
using OptimizationLBFGSB
using Enzyme

# Some stuff to create a posterior object
post # of type Comrade.Posterior

xopt, sol = comrade_opt(post, LBFGSB())
```