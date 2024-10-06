```@meta
CurrentModule = Comrade
```

# Introduction

Comrade is a Bayesian differentiable modular modeling framework for use with very long baseline interferometry.
The goal is to allow the user to easily combine and modify a set of primitive models
to construct complicated source structures. The benefit of this approach is that it is straightforward to construct different source models out of these primitives. Namely, an end-user does
not have to create a separate source "model" every time they
change the model specification. Additionally, most models currently implemented are differentiable with at Enzyme. This allows for gradient accelerated optimization and sampling (e.g., HMC) to be used with little
effort by the end user. 

!!! warning
    As of 0.11 Comrade will only support AD with Enzyme. We have removed support for Zygote and ForwardDiff
    due to performance issues.


## Tutorials

Our tutorial section currently has a large number of examples. The simplest example is fitting simple geometric models to the 2017 M87 data and is detailed in the [Geometric Modeling of EHT Data](@ref) tutorial. We also include "non-parametric" modeling or imaging examples in [Imaging a Black Hole using only Closure Quantities](@ref), and [Stokes I Simultaneous Image and Instrument Modeling](@ref). There is also an introduction to hybrid geometric and image modeling in [Hybrid Imaging of a Black Hole](@ref), which combines physically motivated geometric modeling with the flexibility of image-based models. Finally, we
provide a tutorial on how to use `Comrade` to model polarized data including simultaneously solving for 
the image and instrumental effects like gain ratios and leakage terms in [Polarized Image and Instrumental Modeling](@ref).



## Contributing

This repository has recently moved to [ColPrac](https://github.com/SciML/ColPrac). If you would like to contribute please feel free to open a issue or pull-request.

## Requirements

The minimum Julia version we require is 1.10.
