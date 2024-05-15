# NestedSamples Extension

`Comrade` intefaces to the [`NestedSamplers.jl`](https://github.com/TuringLang/NestedSamplers.jl) package.

We follow `NestedSamplers` interface closely. The 
difference is that instead of creating a `NestedModel`, we 
 pass a `Comrade.VLBIPosterior` object as our model.
Internally, `Comrade` defines the prior transform and extracts the log-likelihood function.

For more information about `NestedSamplers.jl` please see its [docs](https://github.com/TuringLang/NestedSamplers.jl).

## Example

```julia
using Comrade
using NestedSamplers

# Some stuff to create a posterior object
post # of type Comrade.Posterior

# Create sampler using 1000 live points
smplr = Nested(dimension(ascube(post)), 1000)

samples = sample(Random.default_rng(), post, smplr; d_logz=1.0)

# Optionally resample the chain to create an equal weighted output
using StatsBase
equal_weight_chain = Comrade.resample_equal(samples, 10_000)
```