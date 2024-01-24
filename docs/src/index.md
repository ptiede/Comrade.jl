```@meta
CurrentModule = Comrade
```

# Comrade

Comrade is a Bayesian differentiable modular modeling framework for use with very long baseline interferometry.
The goal is to allow the user to easily combine and modify a set of primitive models
to construct complicated source structures. The benefit of this approach is that it is straightforward to construct different source models out of these primitives. Namely, an end-user does
not have to create a separate source "model" every time they
change the model specification. Additionally, most models currently implemented are differentiable with at Zygote and sometimes ForwardDiff[^2]. This allows for gradient accelerated optimization and sampling (e.g., HMC) to be used with little
effort by the end user.
To sample from the posterior, we provide a somewhat barebones interface since, most of the time, and we don't require the additional features offered by most PPLs. Additionally, the overhead introduced by PPLs tends to be rather large. In the future, we may revisit this as
Julia's PPL ecosystem matures.

!!! note
    The primitives the Comrade defines, however, would allow for it to be easily included in PPLs like [`Turing`](https://github.com/TuringLang/Turing.jl).


## Tutorials


## Contributing

This repository has recently moved to [ColPrac](https://github.com/SciML/ColPrac). If you would like to contribute please feel free to open a issue or pull-request.

[^2]: As of 0.9 Comrade switched to using full covariance closures. As a result this requires a sparse cholesky solve in the likelihood evaluation which requires 
a Dual number overload. As a result we recommend using Zygote which does work and
often is similarly performant (reverse 3-6x slower compared to the forward pass).

## Requirements

The minimum Julia version we require is 1.9. In the future we may increase this as Julia advances.


```@contents
Pages = [
    "index.md",
    "vlbi_imaging_problem.md",
    "conventions.md",
    "Tutorials",
    "Libraries",
    "interface.md",
    "base_api.md",
    "api.md"
]
```

## References
