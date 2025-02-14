# AdvancedHMC Extension

The first choice when sampling from the model/image posterior,  is [`AdvancedHMC`](https://github.com/TuringLang/AdvancedHMC.jl), which uses Hamiltonian Monte Carlo to sample from the posterior. Specifically, we usually use the `NUTS` algorithm.

The interface to `AdvancedHMC` is very powerful and general. To simplify 
the procedure for `Stoked` users, we have provided a thin interface. 
A user needs to specify a `sampler` and then call 
the `sample` function.

To sample a user can use follow the standard `AdvancedHMC` interface, e.g.,

```julia
chain = sample(post, NUTS(0.8), 10_000; n_adapts=5_000)
```

!!! warning
    To use HMC the `VLBIPosterior` must be created with a specific `admode` specified.
    The `admode` can be a union of `Nothing` and `<:EnzymeCore.Mode` types. We recommend
    using `Enzyme.set_runtime_activity(Enzyme.Reverse)`


In addition our sample call has a few additional keyword arguments:

 - `saveto = MemoryStore()`: Specifies how to store the samples. The default is `MemoryStore` which stores the samples directly in RAM. For large models this is not a good idea. To save samples periodically to disk use [`DiskStore`](@ref), and then load the results with `load_samples`.

Note that like most `AbstractMCMC` samplers the initial location can be specified with the `initial_params` argument.


## Example 

```julia
using Stoked
using AdvancedHMC
using Enzyme

# Some stuff to create a posterior object
post # of type Stoked.Posterior

out = sample(post, NUTS(0.9), 2_000; n_adapts=1_000, saveto=DiskStore())
chain = load_samples(out)
```
