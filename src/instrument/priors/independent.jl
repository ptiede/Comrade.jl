export IIDSitePrior

abstract type AbstractSitePrior end

segmentation(d::AbstractSitePrior) = getfield(d, :seg)

# ----- the site-prior interface ------------------------------------------------
#
# A site prior describes the prior of one site's instrument parameters over its observed
# (time, frequency) points. `ObservedArrayPrior` routes on `structure`:
#
#   * `ProductStructure()` priors are order-free products over points; when *every* site
#     in an `ArrayPrior` is a product the fast vectorized product-distribution path is
#     used (`build_dist`).
#   * `ChainedStructure()` priors are correlated across their points (e.g. in time).
#     When any site is chained, every site prior — product ones included — contributes
#     an ordinary `Distributions.Distribution` per (site, frequency) chain through
#     `chaindist`; the generic assembly (chained.jl) owns the index maps between the
#     flat storage and the per-chain vectors, the exact conditioning, the merged
#     per-site hyperpriors, and the composed transport nodes.
#
# Extension contract for a new site prior type:
#   required: `structure`, `segmentation`, `chaindist`
#   optional: `hyperprior` (fitted per-site hyperparameters; default none),
#             `anchor_points` (extra exact-conditioning points; default none),
#             `freqseg` (frequency segmentation; default per-channel `SpectralWindow`)
#
# The distribution returned by `chaindist` needs the ordinary `Distributions` API
# (`logpdf`, `rand!`); `Comrade.condition` if it will be combined with reference fixing
# or anchors; and a ProbabilityTransports `transport_node` for the latent transports.
# All three come for free for products of univariate distributions and for
# `GaussMarkovChain`, so wrapping an existing distribution is usually a one-liner.
#
# A future frequency-correlated prior would extend `chaindist` to receive the site's
# whole tuple of chains (with channel centers) and return a single distribution owning
# all of them; the assembly layer is agnostic to how points group into chains.

struct ProductStructure end
struct ChainedStructure end

"""
    structure(sp::AbstractSitePrior) -> ProductStructure() | ChainedStructure()

Declares whether the site prior factorizes as an order-free product over its points
(`ProductStructure`) or is correlated across them (`ChainedStructure`, e.g. priors
correlated in time). Required for every site prior type.
"""
structure(sp::AbstractSitePrior) = throw(
    ArgumentError(
        "$(typeof(sp)) does not implement the Comrade site-prior interface: define " *
            "`Comrade.structure(::$(nameof(typeof(sp)))) = Comrade.ProductStructure()` (or " *
            "`ChainedStructure()`), plus `Comrade.chaindist` if it can appear alongside " *
            "chained site priors."
    )
)

"""
    chaindist(sp::AbstractSitePrior, ts::AbstractVector) -> Distributions.Distribution

Build the distribution of one of the site's chains: the site's parameter values at the
chain times `ts` (in hours, ascending). Sites have one chain per observed frequency
channel; the returned distribution must have `length(d) == length(ts)`. Called once at
`set_array` time.

This is the whole extension surface for a chained site prior: the returned distribution
only needs the ordinary `Distributions` API, plus [`condition`](@ref Comrade.condition)
when used with reference fixing/anchors and a ProbabilityTransports `transport_node` for
the latent transports (all already available for products of univariate distributions
and [`GaussMarkovChain`](@ref Comrade.GaussMarkovChain)).

If the site prior has fitted hyperparameters ([`hyperprior`](@ref Comrade.hyperprior)),
the returned distribution is a *template* holding the hyperparameter placeholders; the
sampled values are threaded through evaluation by the assembly layer (see
`GaussMarkovChain` for the reference implementation).
"""
chaindist(sp::AbstractSitePrior, ts) = throw(
    ArgumentError(
        "$(typeof(sp)) cannot be combined with chained site priors: define " *
            "`Comrade.chaindist(::$(nameof(typeof(sp))), ts)` returning a " *
            "`Distributions.Distribution` over the chain values (see the site-prior " *
            "interface notes in independent.jl)."
    )
)

"""
    hyperprior(sp::AbstractSitePrior) -> NamedTuple

The priors of the site prior's fitted hyperparameters (empty for fully specified
priors). Each site with a non-empty hyperprior gets its own copy, nested under the site
name in the hierarchical sample.
"""
hyperprior(::AbstractSitePrior) = (;)

"""
    anchor_points(sp::AbstractSitePrior, chain::SiteFreqChain) -> (inds, vals)

Extra flat indices of `chain` that the prior exactly conditions on the paired values, on
top of any reference-antenna fixing (e.g. anchored Gauss-Markov chains pin each chain's
first point to zero). Defaults to none.
"""
anchor_points(::AbstractSitePrior, ::SiteFreqChain) = (Int[], Float64[])

"""
    IIDSitePrior(seg::Segmentation, dist; freqseg = SpectralWindow())

Create a site prior that is independent and identically distributed (IID) across all times
and frequencies. The `seg` argument is a segmentation object that defines how fine the time
segmentation is. The `dist` argument is the distribution of the site prior. The `freqseg`
keyword sets the frequency segmentation: the default [`SpectralWindow`](@ref
Comrade.SpectralWindow) gives one independent parameter per observed frequency channel,
while [`FullBand`](@ref Comrade.FullBand) shares a single parameter across the band.

## Example

```julia
A = IIDSitePrior(ScanSeg(), VLBIGaussian(0, 1))
```

creates a site prior that is constant across scans and each scan has a unit Normal prior.

"""
struct IIDSitePrior{S <: Segmentation, F <: FrequencySegmentation, D} <: AbstractSitePrior
    seg::S
    freqseg::F
    dist::D
end

IIDSitePrior(seg::Segmentation, dist; freqseg::FrequencySegmentation = SpectralWindow()) =
    IIDSitePrior(seg, freqseg, dist)

structure(::IIDSitePrior) = ProductStructure()
freqseg(sp::IIDSitePrior) = sp.freqseg

# Deliberately `Product` (not `product_distribution`): packages may fold homogeneous
# vectors into a fused multivariate distribution (e.g. VLBIImagePriors' affine
# pushforwards) whose flat transport is the *centered* identity with logpdf-enforced
# support. The chained assembly needs the per-point bounded transforms and exact
# per-point conditioning that the plain product keeps.
chaindist(sp::IIDSitePrior, ts) = Dists.Product(fill(sp.dist, length(ts)))

"""
    freqseg(sp::AbstractSitePrior) -> FrequencySegmentation

The frequency segmentation of the site prior: how its parameters segment over the
observation's frequency channels. Defaults to [`SpectralWindow`](@ref Comrade.SpectralWindow)
(one independent parameter per channel).
"""
freqseg(::AbstractSitePrior) = SpectralWindow()
