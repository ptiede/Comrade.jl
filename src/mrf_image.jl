export apply_fluctuations, corr_image_prior


"""
    apply_fluctuations(f, mimg::IntensityMap, δ::AbstractArray)

Apply multiplicative fluctuations to an image `mimg` with fluctuations `δ`.
The function `f` is applied to the fluctuations and then the the transfored δ are multiplicatively applied
to the image.
"""
@inline function apply_fluctuations(f, mimg::IntensityMap, δ::AbstractArray)
    return IntensityMap(_apply_fluctuations(f, baseimage(mimg), δ), axisdims(mimg))
end

@inline function apply_fluctuations(f, m::AbstractModel, g::AbstractRectiGrid, δ::AbstractArray)
    return apply_fluctuations(f, intensitymap(m, g), δ)
end

@inline function apply_fluctuations(t::VLBIImagePriors.LogRatioTransform, m::AbstractModel, g::AbstractRectiGrid, δ::AbstractArray)
    mimg = baseimage(intensitymap(m, g))
    return apply_fluctuations(t, IntensityMap(mimg ./ sum(mimg), g), δ)
end


@inline function apply_fluctuations(mimg::IntensityMap, δ::AbstractArray)
    return apply_fluctuations(identity, mimg, δ)
end

@inline function _apply_fluctuations(f, mimg::AbstractArray, δ::AbstractArray)
    return mimg .* f.(δ)
end

@noinline _checknorm(m::AbstractArray) = isapprox(sum(m), 1, atol = 1.0e-6)
EnzymeRules.inactive(::typeof(_checknorm), args...) = nothing

function _fastsum(x)
    tot = zero(eltype(x))
    @simd for i in eachindex(x)
        tot += x[i]
    end
    return tot
end


function _apply_fluctuations(t::VLBIImagePriors.LogRatioTransform, mimg::AbstractArray, δ::AbstractArray)
    @argcheck _checknorm(mimg) "Mean image must have unit flux when using log-ratio transformations in apply_fluctuations while it seems to be $(sum(mimg))"
    r = to_simplex(t, baseimage(δ))
    r .= r .* baseimage(mimg)
    r .= r ./ _fastsum(r)
    return r
end


"""
    corr_image_prior(grid::AbstractRectiGrid, corr_length::Real; base=GMRF, order=1, lower=1.0, upper=2*max(size(grid)...))
    corr_image_prior(grid::AbstractRectiGrid, obs::EHTObservationTable; base=GMRF, order=1, lower=1.0, upper=2*max(size(grid)...))

Construct a correlated image prior, for the image with grid `grid`, and using the observation `dvis`. 
The correlation will be a Markov Random Field (MRF) of order `order`, with the base distribution `base`. For 
`base` you can choose any of the Markov random fields defined in `VLBIImagePriors`, the default is `GMRF` which
is a Gaussian MRF. 

As part the prior will be a hierarchical prior with the correlation length as a hyperparameter. By default the correlation
parameter uses a first order inverse gamma distribution for its prior. The `frac_below_beam` parameter is the fraction of the
correlation prior mass that is below the beam size of the observation `dvis`. The `lower` and `upper` parameters are the lower
and upper bounds of the correlation length, we don't let the correlation length to be too small or large for numerical reasons.

## Arguments
 - `grid::AbstractRectiGrid`: The grid of the image to be reconstructed.
 - `corr_length`: The correlation length of the MRF. If this is an `EHTObservationTable` then the corr_length 
    will be the approximate beam size of the observation. 

## Keyword Arguments
 - `base`: The base distribution of the MRF. Options include `GMRF`, `EMRF`, and `CMRF`
 - `order`: The order of the MRF. Default is first order
 - `frac_below_beam`: The fraction of the correlation prior mass that is below the beam size of the observation `dvis`. 
    the default is `0.01` which means only 1% of the log-image correlation length is below the beam size.
 - `lower`: The lower bound of the correlation length. Default is `1.0`
 - `upper`: The upper bound of the correlation length. Default is `2*max(size(grid)...)`


!!! warn
    An order > 2 will be slow since we switch to a sparse matrix representation of the MRF.
"""
function corr_image_prior(
        grid::AbstractRectiGrid, corr_length::Real;
        base = GMRF, order = 1,
        frac_below_beam = 0.01,
        lower = 1.0, upper = 2 * max(size(grid)...)
    )
    rat = corr_length / step(grid.X)
    cmarkov = ConditionalMarkov(base, grid; order = order)
    dρ = Dists.truncated(Dists.InverseGamma(1.0, -log(frac_below_beam) * rat); lower, upper)
    cprior = HierarchicalPrior(cmarkov, dρ)
    return cprior
end

corr_image_prior(grid::AbstractRectiGrid, obs::EHTObservationTable; kwargs...) = corr_image_prior(grid, beamsize(obs); kwargs...)
