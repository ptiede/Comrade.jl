```@meta
CurrentModule = Comrade
```

# Comrade

Comrade is a **differentiable** modular modeling framework for use with very long baseline interferometry.
The goal is to allow the user to easily combine and modify a set of *primitive* models
to construct complicated source structures. The benefit of this approach is that is straightforward to construct different source models out of these primitives. Namely, a end-user does
not have to create a separate source "model" everytime they 
change the model specification. Additionally, most models currently implemented are differentiable with at least `ForwardDiff`. This allows for gradient accelerated optimization, and sampling (e.g. HMC) to be used with little
effort by the end user. 


`Comrade` does not currently have a native optimization or 
sampling interface. The reasoning for this is that different 
problems are amenable to different optimizers. Rather than 
including all optimizers in Comrade, expanding the number of
dependencies, Comrade tries to make moving from an image model
to objective function easy. As an example of this we 
To use perform inferences on data you can then hook into the vast array of different 
modeling and optimization packages in Julia. There are some small examples packages
defining these interface such as [ComradeSoss.jl](https://github.com/ptiede/ComradeSoss.jl) 
which combines Comrade with `Soss` a probabilistic programming language. Other interfaces
to e.g. [Turing](https://turing.ml/stable/), [BAT](https://github.com/bat/BAT.jl) are 
planned.

## Requirements

The minimum Julia version we require is 1.6, which is the current LTS release. In the 
future we may increase this as Julia advances.

```@contents
Pages = [
    "index.md",
    "vlbi_imaging_problem.md",
    "example.md",
    "interface.md",
    "api.md"
]
```
