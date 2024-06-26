export apply_fluctuations


"""
    apply_fluctuations(f, mimg::IntensityMap, δ::AbstractArray)

Apply multiplicative fluctuations to an image `mimg` with fluctuations `δ`.
The function `f` is applied to the fluctuations and then the the transfored δ are multiplicatively applied
to the image.
"""
function apply_fluctuations(f, mimg::IntensityMap, δ::AbstractArray)
    return IntensityMap(_apply_fluctuations(f, baseimage(mimg), δ), axisdims(mimg))
end

function apply_fluctuations(f, m::AbstractModel, g::AbstractRectiGrid, δ::AbstractArray)
    return apply_fluctuations(f, intensitymap(m, g), δ)
end

function apply_fluctuations(t::VLBIImagePriors.LogRatioTransform, m::AbstractModel, g::AbstractRectiGrid, δ::AbstractArray)
    # Hack to prevent Zygote from trying to AD through IntensityMap constants
    mimg = baseimage(intensitymap(m, g))
    return apply_fluctuations(t, IntensityMap(mimg./sum(mimg), g), δ)
end



function apply_fluctuations(mimg::IntensityMap, δ::AbstractArray)
    return apply_fluctuations(identity, mimg, δ)
end

function _apply_fluctuations(f, mimg::AbstractArray, δ::AbstractArray)
    return mimg.*f.(δ)
end

function _apply_fluctuations(t::VLBIImagePriors.LogRatioTransform, mimg::AbstractArray, δ::AbstractArray)
    @argcheck isapprox(sum(mimg), 1, atol=1e-6) "Mean image must have unit flux when using log-ratio transformations in apply_fluctuations"
    r = baseimage(mimg).*to_simplex(t, δ)
    return r./sum(r)
end
