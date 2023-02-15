```@meta
CurrentModule = Comrade
```

# Comrade

Comrade is a Bayesian differentiable modular modeling framework for use with very long baseline interferometry.
The goal is to allow the user to easily combine and modify a set of primitive models
to construct complicated source structures. The benefit of this approach is that it is straightforward to construct different source models out of these primitives. Namely, an end-user does
not have to create a separate source "model" every time they
change the model specification. Additionally, most models currently implemented are differentiable with at least ForwardDiff and Zygote. This allows for gradient accelerated optimization and sampling (e.g., HMC) to be used with little
effort by the end user.
To sample from the posterior, we provide a somewhat barebones interface since, most of the time, and we don't require the additional features offered by most PPLs. Additionally, the overhead introduced by PPLs tends to be rather large. In the future, we may revisit this as
Julia's PPL ecosystem matures.

!!! note
    The primitives the Comrade defines, however, would allow for it to be easily included in PPLs like [Turing](https://github.com/TuringLang/Turing.jl.



Our tutorial section currently has a large number of examples. The simplest example is fitting simple geometric models to the 2017 M87 data and is detailed in the [Making an Image of a Black Hole](@ref) tutorial. We also include "non-parametric" modeling or imaging examples in [Imaging a Black Hole using only Closure Quantities](@ref), and
[Stokes I Simultaneous Image and Instrument Modeling](@ref). There is also an introduction to hybrid geometric and image modeling in [Hybrid Imaging of a Black Hole](@ref), which combines physically motivated geometric modeling with the flexibility of image-based models.


As of 0.7, Comrade also can simultaneously reconstruct polarized image models and instrument corruptions through the RIME[^1] formalism. A short example explaining
these features can be found in [Polarized Image and Instrumental Modeling](@ref).
## Contributing

This repository has recently moved to [ColPrac](https://github.com/SciML/ColPrac). If you would like to contribute please feel free to open a issue or pull-request.


## Requirements

The minimum Julia version we require is 1.7. In the future we may increase this as Julia advances.


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
[^1]: Hamaker J.P and Bregman J.D. and Sault R.J. Understanding radio polarimetry. I. Mathematical foundations [ADS](https://ui.adsabs.harvard.edu/abs/1996A&AS..117..137H). 

