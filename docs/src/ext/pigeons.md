# Pigeons Extension

`Comrade` intefaces to the [`Pigeons.jl`](https://github.com/Julia-Tempering/Pigeons.jl) package.

The interface to Pigeons follows the same one as in the Pigeons docs. Since `Comrade` implements the
`LogDensityProblems` interface it can be used as is. Note that the user must transform the target
to either Rⁿ or the unity hypercube using `asflat` or `ascube` respectively, before passing it to `Pigeons`.

Additionally, `pigeons` will return the standard output. If you want to directly access the samples you can use
the `sample_array` function and pass the transformed target and the samples object to convert them to the usual
parameter space.

For more information about `Pigeons.jl` please see its [docs](https://pigeons.run/dev/).

## Example

```julia
using Comrade
using Pigeons

# Some stuff to create a posterior object
post # of type Comrade.Posterior

# Create sampler using 1000 live points
samples = pigeons(target=ascube(post), explorer=SliceSampler(), record=[traces])
# Transform the samples to the parameter space
chain = sample_array(ascube(post), samples)
```