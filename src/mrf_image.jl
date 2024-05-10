
"""
    norm_image_noise(mean, rf, g)

Creates a image model with some mean image `mean` and a random field `rf` whose
noise is added multiplicatively to the mean image. The grid `g` is used to
specify the domain of the image.
"""
function norm_image_noise(mean::AbstractArray, rf::AbstractArray, g::AbstractRectiGrid)
    addn =  mean./sum(mean) .* rf./sum(rf)
    return IntensityMap(addn./sum(addn), g)
end

function norm_image_noise(mean::AbstractModel, rf::AbstractArray, g::AbstractRectiGrid)
    mimg = intensitymap(mean, g)
    return image_noise(baseimage(mimg), rf, g)
end
