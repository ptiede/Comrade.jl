# ComradeAHMC

The first choice when sampling from the model/image posterior,  is [`AdvancedHMC`](https://github.com/TuringLang/AdvancedHMC.jl)
), which uses Hamiltonian Monte Carlo to sample from the posterior. Specifically, we usually use the `NUTS` algorithm. 

The interface to `AdvancedHMC` is very powerful and general. To simplify 
the procedure for `Comrade` users, we have provided a thin interface. 
A user needs to specify a `sampler` and then call 
the `sample` function.

For `AdvancedHMC`, the user can create the sampler by calling the [`AHMC`](@ref) function. This only has one mandatory argument, the `metric` the sampler uses. There are currently two options:

    - `DiagEuclideanMetric` which uses a diagonal metric for covariance adaptation
    - `DenseEuclideanMetric` which uses a dense or full rank metric for covariance adaptation

We recommend that a user starts with `DiagEuclideanMetric` since the dense metric typically requires many more samples to tune correctly. 
The other options for `AHMC` (sans `autodiff`) specify which version of HMC to use. Our default options match the choices made by the [Stan](https://mc-stan.org/) programming language. The final option to consider is the `autodiff` optional argument. This specifies which auto differentiation package to use. Currently
`Val(:Zygote)` is the recommended default for all models. If you model doesn't
work with Zygote please file an issue. Eventually we will move entirely to Enzyme.

## Example 

```julia
using Comrade
using ComradeAHMC

# Some stuff to create a posterior object
post # of type Comrade.Posterior

metric = DiagEuclideanMetric(dimension(post))
smplr = AHMC(metric=metric, autodiff=Val(:Zygote))

samples, stats = sample(post, smplr, 2_000; nadapts=1_000)
```

## API

```@meta
CurrentModule = ComradeAHMC
```

```@autodocs
Modules = [ComradeAHMC]
```

