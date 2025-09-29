export apply_fluctuations, apply_fluctuations!, UnitFluxMap, corr_image_prior


"""
    apply_fluctuations(f, mimg::IntensityMap, δ::AbstractArray)

Apply multiplicative fluctuations to an image `mimg` with fluctuations `δ`.
The function `f` is applied to the fluctuations and then the the transfored δ are multiplicatively applied
to the image.
"""
@inline function apply_fluctuations(f, mimg::IntensityMap, δ::AbstractArray)
    out = similar(mimg)
    apply_fluctuations!(f, out, mimg, δ)
    return out
end

@inline function apply_fluctuations(f, m::AbstractModel, g::AbstractRectiGrid, δ::AbstractArray)
    return apply_fluctuations(f, intensitymap(m, g), δ)
end

@inline function apply_fluctuations(t::VLBIImagePriors.LogRatioTransform, m::AbstractModel, g::AbstractRectiGrid, δ::AbstractArray)
    mimg = baseimage(intensitymap(m, g))
    return apply_fluctuations(t, IntensityMap(mimg ./ _fastsum(mimg), g), δ)
end


@inline function apply_fluctuations(mimg::IntensityMap, δ::AbstractArray)
    return apply_fluctuations(identity, mimg, δ)
end

@inline function apply_fluctuations!(f, out::IntensityMap, mimg::IntensityMap, δ::AbstractArray)
    bout = baseimage(out)
    bmimg = baseimage(mimg)
    bout .= bmimg .* f.(δ)
    return nothing
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


"""
    UnitFluxMap(f)

A transformation that broadcasts a function `f` over an array, and the normalizes the resulting array
to have unit flux. This is used with the function [`apply_fluctuations`](@ref) while 
imaging.
"""
struct UnitFluxMap{F}
    f::F
end


"""
    apply_fluctuations!(t::UnitFluxMap{F}, out::IntensityMap, mimg::IntensityMap, δ::AbstractArray)

Apply multiplicative fluctuations to an image `mimg` with fluctuations `δ` using the transformation `t`,
by broadcasting the function `F` over the array `δ` and then normalizing the resulting array to have unit flux.

"""
function apply_fluctuations!(t::UnitFluxMap, out::IntensityMap, mimg::IntensityMap, δ::AbstractArray)
    @argcheck _checknorm(mimg) "Mean image must have unit flux when using unit flux transformations in apply_fluctuations while it seems to be $(sum(mimg))"
    f = t.f
    bout = baseimage(out)
    @inbounds for i in eachindex(bout, δ)
        bout[i] = f(δ[i])
    end
    # bout .= f.(baseimage(δ))
    fd = _fastsum(bout)
    bmimg = baseimage(mimg)

    @inbounds for i in eachindex(bout, bmimg)
        bout[i] *= bmimg[i] / fd
    end

    fi = _fastsum(bout)
    bout .*= inv(fi)

    return nothing
end


function apply_fluctuations!(t::VLBIImagePriors.LogRatioTransform, out::IntensityMap, mimg::IntensityMap, δ::AbstractArray)
    @argcheck _checknorm(mimg) "Mean image must have unit flux when using log-ratio transformations in apply_fluctuations while it seems to be $(sum(mimg))"
    bout = baseimage(out)
    to_simplex!(t, bout, baseimage(δ))

    bmimg = baseimage(mimg)
    for i in eachindex(bout, bmimg)
        bout[i] *= bmimg[i]
    end
    fi = _fastsum(bout)

    bout .*= inv(fi)
    return nothing
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
