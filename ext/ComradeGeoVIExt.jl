module ComradeGeoVIExt

using Comrade
using GeoVI
import VLBILikelihoods
import GeoVI: logdensity, normalized_residual, transformation,
    leftsqrtmetric, rightsqrtmetric, fishermetric, compose

# ----- the existing Comrade/VLBILikelihoods visibility likelihood AS a GeoVI likelihood ----
#
# We deliberately do NOT introduce a bespoke `<: GeoVI.AbstractLikelihood` wrapper. GeoVI's
# `compose` accepts any object that implements the observation-space likelihood interface, so
# we implement that interface directly on `VLBILikelihoods.ComplexVisLikelihood` — the same
# diagonal complex-Gaussian likelihood Comrade already builds in `makelikelihood`. `d.μ` holds
# the data visibilities and `d.Σ` the per-visibility variance (σ²); the prediction `ŷ` is the
# model visibilities. The log density reuses VLBILikelihoods' own `unnormed_logpdf` (NaN-safe
# for flagged data); the Fisher metric of the diagonal Gaussian is elementwise `1/Σ`, so its
# square-root actions are elementwise `1/√Σ`. (The standard-normal prior is supplied by GeoVI.)

const CVL = VLBILikelihoods.ComplexVisLikelihood

logdensity(d::CVL, ŷ)              = VLBILikelihoods.unnormed_logpdf(d, ŷ)
normalized_residual(d::CVL, ŷ)     = (d.μ .- ŷ) ./ sqrt.(d.Σ)
transformation(d::CVL, ŷ)          = ŷ ./ sqrt.(d.Σ)
leftsqrtmetric(d::CVL, ŷ, η)       = η ./ sqrt.(d.Σ)
rightsqrtmetric(d::CVL, ŷ, v)      = v ./ sqrt.(d.Σ)
fishermetric(d::CVL, ŷ, v)         = v ./ d.Σ

# ----- assemble the composed GeoVI likelihood from a standardized posterior ---------------

function Comrade.geovi_likelihood(tpost::Comrade.TransformedVLBIPosterior;
        adtype = GeoVI.AutoEnzyme())
    post = tpost.lpost
    data = dataproducts(post)
    length(data) == 1 || throw(ArgumentError(
        "geovi_likelihood currently supports a single visibility data product; got $(length(data))"))
    dvis = first(data)
    # the same construction Comrade's `makelikelihood` uses — reuse the data-likelihood type
    lh = VLBILikelihoods.ComplexVisLikelihood(measurement(dvis), noise(dvis) .^ 2)
    forward = Base.Fix1(forwardmodel, tpost)
    return compose(lh, forward; adtype = adtype)
end

function forwardmodel(tpost, ξ)
    # the same construction Comrade's `makelikelihood` uses — reuse the forward model
    return last(forward_model(tpost.lpost, transform(tpost, ξ)))
end

# A default starting latent for the `geovi_likelihood` composed likelihood, so GeoVI's
# `init`/`fit` can be called without an explicit `ξ0`. Size = `dimension(tpost)`; the array
# TYPE (host vs Reactant device) is taken from the data-likelihood's variance array `Σ`, so the
# returned zero latent matches the problem's backend. Dispatch keys on the forward being
# Comrade's `forwardmodel`, so this only fires for `geovi_likelihood`-built likelihoods.
function GeoVI.default_latent(
        lh::GeoVI.ComposedLikelihood{L, F, LIN, AD, AO},
    ) where {L <: CVL, F <: Base.Fix1{typeof(forwardmodel)}, LIN, AD, AO}
    tpost = lh.forward.x
    Σ = lh.likelihood.Σ
    return fill!(similar(Σ, real(eltype(Σ)), Comrade.dimension(tpost)), 0)
end

end # module
