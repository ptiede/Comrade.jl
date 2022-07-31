# ComradeAHMC

To sample from the posterior our first choice is typically [`AdvancedHMC`](https://github.com/TuringLang/AdvancedHMC.jl)
) which uses Hamiltonian Monte Carlo to sample from the posterior. Specifically we usually use the `NUTS` algorithm. 

The interface to `AdvancedHMC` is very powerful and general. To simplify 
the procedure for `Comrade` users we have provided a thin interface. 
Like all `Comrade` a user just needs to specify a `sampler` then call 
the `sample` function.

For `AdvancedHMC` the user can create the sampler by calling the [`AHMC`](@ref) function. This only has one mandatory argument the `metric` the sampler uses. There are typically two options:

    - `DiagEuclideanMetric` which uses a diagonal metric for covariance adaptation
    - `DenseEuclideanMetric` which uses a dense or full rank metric for covariance adaptation

We typically recommend a user to start with `DiagEuclideanMetric` since the dense metric typically requires a lot more samples to properly tune. 
The other options for `AHMC` (sans `autodiff`) specify exactly which version of HMC to use. Our default options match the choices made by the [Stan](https://mc-stan.org/) programming language. The final option to consider is the `autodiff` optional argument. This specifies which autodifferentiation package to use. For geometric modeling we recommend the `AD.ForwardDiffBackend()`, while for Bayesian Imaging `AD.ZygoteBackend()`. 

## Example 

```julia
using Comrade
using ComradeAHMC

# Some stuff to create a posterior object
post # of type Comrade.Posterior

metric = DiagEuclideanMetric(dimension(post))
smplr = AHMC(metric=metric, autodiff=AD.ForwardDiffBackend{10}())

chain, stats = sample(post, smplr, 2_000; nadapts=1_000)
```

## API

```@meta
CurrentModule = ComradeAHMC
```

```@autodocs
Modules = [ComradeAHMC]
```

