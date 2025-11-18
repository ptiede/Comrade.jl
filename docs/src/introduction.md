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

Our tutorial section currently has a large number of examples. Including:
 - [Geometric Modeling of EHT Data](@ref)
 - [Imaging a Black Hole using only Closure Quantities](@ref)
 - [Stokes I Simultaneous Image and Instrument Modeling](@ref)
 - [Hierarchicial Interferometric Bayesian Imaging (HIBI)](@ref)
 - [Modeling the Power Spectrum of an AGN with Markov Random Field Expansion](@ref)
 - [Polarized Image and Instrumental Modeling](@ref)

Note that these tutorials should be treated as starting points for your own modeling efforts! 
`Comrade` is highly flexible and so you should be able to combine different aspects of these 
tutorials to create your own models. We are currently working on a high-level imaging 
code that uses `Comrade` as a backend, but in the mean time these tutorials should help you get 
on your own VLBI imaging projects.


## Contributing

This repository has recently moved to [ColPrac](https://github.com/SciML/ColPrac). If you would like to contribute please feel free to open a issue or pull-request.

## Requirements

The minimum Julia version we require is 1.10.
