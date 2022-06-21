```@meta
CurrentModule = Comrade
```

# Comrade

Comrade is a **Bayesian differentiable** modular modeling framework for use with very long baseline interferometry.
The goal is to allow the user to easily combine and modify a set of *primitive* models
to construct complicated source structures. The benefit of this approach is that is straightforward to construct different source models out of these primitives. Namely, a end-user does
not have to create a separate source "model" every time they
change the model specification. Additionally, most models currently implemented are differentiable with at least `ForwardDiff`. This allows for gradient accelerated optimization, and sampling (e.g. HMC) to be used with little
effort by the end user.

Currently there are two main ways to construct and sample from the posterior. The first is the simple but somewhat more limited
native interface. To see how this works see the [Making an Image of a Black Hole](@ref) tutorial. The other method is [ComradeSoss.jl](https://github.com/ptiede/ComradeSoss.jl) which combines Comrade with `Soss` a probabilistic programming language. This allows for easier composition of models, and provides a more complete Bayesian workflow, including the ability
to sample from the posterior predictive distributions. Other interfaces to e.g. [Turing](https://turing.ml/stable/), [BAT](https://github.com/bat/BAT.jl) are planned.

## Contributing

This repository has recently moved to [ColPrac](https://github.com/SciML/ColPrac). If you would like to contribute please feel free to open a issue or pull-request.


## Requirements

The minimum Julia version we require is 1.6, which is the current LTS release. In the future we may increase this as Julia advances.

```@contents
Pages = [
    "index.md",
    "vlbi_imaging_problem.md",
    "Tutorials" => [
    "examples/data.md",
    "examples/black_hole_image.md",
    "examples/nonanalytic.md"
    ]
    "interface.md",
    "base_api.md",
    "api.md"
]
```
