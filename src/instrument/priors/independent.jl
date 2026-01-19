export IIDSitePrior, RWSitePrior

"""
    AbstractSitePrior

An abstract type for site priors. This defines an abstract distribution for the prior of a single
site (antenna). Specific site priors should subtype this abstract type and then implement 
```julia
sitedist(d::AbstractSitePrior, site::Symbol, time, frequency, sitemap::SiteMap)
```
"""
abstract type AbstractSitePrior end

"""
    sitedist(d::AbstractSitePrior, site::Symbol, time, frequency, smap)
Get the distribution for a specific site at a specific time and frequency,
"""
function sitedist end

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

"""
    RWSitePrior(seg::Segmentation, dist0, transition)

Create a site prior that is essentially a random walk from segment to segment.
The `seg` argument is a segmentation object that defines how fine the time
segmentation is. The `dist0` argument is the distribution for the initial scan,
and the `transition` argument is the distribution for the transition between segments.

This means
x0 ~ dist0
xt = x0 + ϵt,  where ϵt ~ transition


## Example

```julia
A = RWSitePrior(ScanSeg(),  Normal(0, 1), Normal(0, 0.1))
```

"""
struct RWSitePrior{S <: Segmentation, D0, DT} <: AbstractSitePrior
    seg::S
    dist0::D0
    trans::DT
end
