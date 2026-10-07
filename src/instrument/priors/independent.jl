export IIDSitePrior

abstract type AbstractSitePrior end

segmentation(d::AbstractSitePrior) = getfield(d, :seg)

# Whether this site prior's observed distribution needs the per-chain machinery of
# gauss_markov.jl rather than a product distribution (see `ObservedArrayPrior`).
# Time-correlated site-prior kinds opt in by overriding this.
needs_chain_machinery(::AbstractSitePrior) = false

"""
    IIDSitePrior(seg::Segmentation, dist; init = nothing)

Create a site prior that is independent and identically distributed (IID) across all times
and frequencies. The `seg` argument is a segmentation object that defines how fine the time
segmentation is. The `dist` argument is the distribution of the site prior.

`init = FixedInit(v)` fixes each site's *first* time stamp (per frequency channel) to `v`
instead of sampling it — the IID analogue of the chain-opening conditioning of
[`GaussMarkovSitePrior`](@ref). Use it when the parameter rides on a separate per-site
offset term: without the pin the offset trades off against a common shift of the site's
values, and for a phase the likelihood's `2π` periodicity turns that redundancy into
spurious modes. A stamp the referencing scheme (`refant`) already fixes stays the
reference's; the two values must agree.

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

# The first-stamp pin of this site prior, or `nothing`. `build_dist` folds it into the
# fixed indices alongside the referencing scheme's.
initprior(::AbstractSitePrior) = nothing
initprior(d::IIDSitePrior) = d.init
