# ComradeAdaptMCMC

Experimental interface to the [`AdaptiveMCMC.jl](https://github.com/mvihola/AdaptiveMCMC.jl) MCMC package. This uses parallel tempering to sample from the posterior. We typically recommend using one of the nested sampling packages.
This interface follows `Comrade`'s usual sampling interface for uniformity.

## Example

```julia
using Comrade
using ComradeAdaptMCMC

# Some stuff to create a posterior object
post # of type Comrade.Posterior


smplr = AdaptMCMC(ntemp=5) # use 5 tempering levels

chain, stats = sample(post, smplr, 500_000, 300_000)
```

## API

```@meta
CurrentModule = ComradeAdaptMCMC
```

```@autodocs
Modules = [ComradeAdaptMCMC]
```