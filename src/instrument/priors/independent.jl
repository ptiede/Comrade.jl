export IIDSitePrior

abstract type AbstractSitePrior end

segmentation(d::AbstractSitePrior) = getfield(d, :seg)

"""
    IIDSitePrior(seg::Segmentation, dist)

Create a site prior that is independent and identically distributed (IID) across all times
and frequencies. The `seg` argument is a segmentation object that defines how fine the time
segmentation is. The `dist` argument is the distribution of the site prior.

## Example

```julia
A = IIDSitePrior(ScanSeg(), Normal(0, 1))
```

creates a site prior that is constant across scans and each scan has a unit Normal prior.

"""
struct IIDSitePrior{S <: Segmentation, D} <: AbstractSitePrior
    seg::S
    dist::D
end
