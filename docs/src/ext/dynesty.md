# Dynesty Extension

`Dynesty` interfaces `Comrade` to the excellent [`dynesty`](https://github.com/joshspeagle/dynesty) package, more specifically
the [Dynesty.jl](https://github.com/ptiede/Dynesty.jl) Julia wrapper.

We follow `Dynesty.jl` interface closely. However, 
instead of having to pass a log-likelihood function and prior transform, we instead just pass a `Comrade.VLBIPosterior` object
and `Comrade` takes care of defining the prior transformation and 
log-likelihood for us.
For more information about `Dynesty.jl`, please see its [docs](https://github.com/ptiede/Dynesty.jl) and docstrings.

## Example

```julia
using Comrade
using Dynesty

# Some stuff to create a posterior object
post # of type Comrade.Posterior

# Create sampler using 1000 live points
smplr = NestedSampler(;nlive=1000)

chain = dysample(post, smplr; dlogz=1.0)

# Optionally resample the chain to create an equal weighted output
using StatsBase
equal_weight_chain = Comrade.resample_equal(samples, 10_000)
```
