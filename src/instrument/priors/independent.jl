export IIDSitePrior

abstract type AbstractSitePrior end

segmentation(d::AbstractSitePrior) = getfield(d, :seg)

# Whether the observed prior is built by gauss_markov.jl rather than as a product
# distribution. Time-correlated site priors return `true`.
needs_chain_machinery(::AbstractSitePrior) = false

"""
    IIDSitePrior(seg::Segmentation, dist; init = nothing)

Create a site prior that is independent and identically distributed (IID) across all times
and frequencies. The `seg` argument is a segmentation object that defines how fine the time
segmentation is. The `dist` argument is the distribution of the site prior.

`init = FixedInit(v)` fixes each site's first time stamp (per frequency channel) to `v`
instead of sampling it, as the same option does for [`GaussMarkovSitePrior`](@ref).
Use it when the parameter is added to a separate per-site offset term: otherwise the offset
and a common shift of the site's values are degenerate, and for a phase the `2π`
periodicity of the likelihood turns that degeneracy into extra modes. A stamp that the
referencing scheme (`refant`) also fixes must get the same value.

## Example

```julia
A = IIDSitePrior(ScanSeg(), VLBIGaussian(0, 1))
```

creates a site prior that is constant across scans and each scan has a unit Normal prior.

"""
struct IIDSitePrior{S <: Segmentation, D, I} <: AbstractSitePrior
    seg::S
    dist::D
    init::I
end

function IIDSitePrior(seg::Segmentation, dist; init = nothing)
    init === nothing || init isa FixedInit || throw(
        ArgumentError(
            "IIDSitePrior only supports `init = FixedInit(v)` (or nothing): there is " *
                "no chain whose start a distributional init could describe"
        )
    )
    return IIDSitePrior(seg, dist, init)
end

# The first-stamp pin of this site prior, or `nothing`.
initprior(::AbstractSitePrior) = nothing
initprior(d::IIDSitePrior) = d.init
