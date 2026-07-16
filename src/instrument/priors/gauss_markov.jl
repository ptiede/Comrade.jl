export GaussMarkovSitePrior

"""
    GaussMarkovSitePrior(seg::TimeSegmentation, process::AbstractGaussMarkovProcess; centered = false, anchored = false)

A site prior that is correlated in time, following the continuous-time Gauss-Markov
`process` (e.g. [`OrnsteinUhlenbeck`](@ref), [`Wiener`](@ref)) sampled at the times
implied by the segmentation `seg`. The log-density factorizes over each site's time chain
using the process's exact discretization, so the cost is `O(n)` in the number of times and
the irregular time sampling of a VLBI observation is handled exactly. When the observation
has multiple frequency channels each (site, channel) pair is an independent chain.

Process fields that are `Distributions.Distribution`s are *fitted* hyperparameters. When
any are present the prior sample is a `NamedTuple` `(params = SiteArray, hyperparams =
NamedTuple)`, mirroring `HierarchicalPrior`; when all fields are numbers the sample is a
plain `SiteArray` like other site priors. Every site gets its *own* hyperparameters: the
process passed to `ArrayPrior` acts as a per-site template, and each site's fitted fields
are nested under the site's name (e.g. `x.instrument.lg.hyperparams.AA.σ`). Per-site
overrides via the usual `ArrayPrior` kwargs change the template (or the prior entirely)
for that site. `caltable`/`plotcaltable` accept the hierarchical sample directly; for
elementwise postprocessing (e.g. `exp.(...)` of log-gains) use the `params` field, i.e.
`exp.(x.instrument.lg.params)`.

By default the prior is *whitened* (non-centered): each free point's latent coordinate is
an iid standard variate that is colored through the chain's exact conditional moments.
This removes the hyperparameter-gain funnel when hyperparameters are fitted and makes the
transport to the Std spaces exact, so both `asflat` and `ascube` work. Set
`centered = true` to instead use the gain values themselves as coordinates (the centered
parameterization), which can mix better when every point is strongly data-constrained;
centered priors support `asflat` only.

This prior is intended for both gain amplitudes and phases. Phases are modeled as
*unwrapped* reals (the likelihood only uses `exp(iθ)`, so wrapping is a display concern).
Reference stations (`refant`) are handled by exact conditioning of the chain on the fixed
values, which works with scattered fixed indices such as those produced by
`SEFDReference`.

Setting `anchored = true` additionally conditions each site's chain to be exactly zero at
its first time, so the prior describes the *evolution* of a quantity relative to its
starting value. This is important when the chain models an unwrapped phase fluctuation
on top of a separate offset term: without anchoring the chain's overall level trades off
against the offset, and the likelihood's `2π` periodicity turns that redundancy into
spurious posterior modes at `2π`-shifted levels. Anchoring removes the level freedom
entirely, so the wrap ambiguity lives only in the (circular) offset while genuine
continuous drift is still expressed by the chain. Processes without a proper initial law
([`Wiener`](@ref)) require `anchored = true` (or a reference scheme fixing each chain's
first point).

Under the hood this simply builds a [`GaussMarkovChain`](@ref Comrade.GaussMarkovChain)
per (site, frequency) chain — see that distribution for the standalone chain semantics.

# Example

```julia
lg ~ ArrayPrior(GaussMarkovSitePrior(IntegSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0))))
gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 3.0, τ = 1.0)); refant = SEFDReference(0.0))
## phase-fluctuation term meant to be combined with a circular per-track offset
dgp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), Wiener(σ = 1.0); anchored = true); refant = SEFDReference(0.0))
```
"""
struct GaussMarkovSitePrior{S <: TimeSegmentation, F <: FrequencySegmentation, P <: AbstractGaussMarkovProcess} <: AbstractSitePrior
    seg::S
    freqseg::F
    process::P
    centered::Bool
    anchored::Bool
end

function GaussMarkovSitePrior(
        seg::TimeSegmentation, process::AbstractGaussMarkovProcess;
        freqseg::FrequencySegmentation = SpectralWindow(), centered = false, anchored = false
    )
    return GaussMarkovSitePrior(seg, freqseg, process, centered, anchored)
end

structure(::GaussMarkovSitePrior) = ChainedStructure()
freqseg(sp::GaussMarkovSitePrior) = sp.freqseg
hyperprior(sp::GaussMarkovSitePrior) = hyperprior(sp.process)

function anchor_points(sp::GaussMarkovSitePrior, c::SiteFreqChain)
    sp.anchored || return (Int[], Float64[])
    return ([first(c.inds)], [0.0])
end

chaindist(sp::GaussMarkovSitePrior, ts::AbstractVector) =
    GaussMarkovChain(sp.process, ts; centered = sp.centered)
