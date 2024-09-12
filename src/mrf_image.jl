export apply_fluctuations


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
    return apply_fluctuations(t, IntensityMap(mimg./sum(mimg), g), δ)
end



@inline function apply_fluctuations(mimg::IntensityMap, δ::AbstractArray)
    return apply_fluctuations(identity, mimg, δ)
end

@inline function _apply_fluctuations(f, mimg::AbstractArray, δ::AbstractArray)
    return mimg.*f.(δ)
end

@noinline _checknorm(m::AbstractArray) = isapprox(sum(m), 1, atol=1e-6)
Enzyme.EnzymeRules.inactive(::typeof(_checknorm), args...) = nothing

function _fastsum(x)
    tot = zero(eltype(x))
    @simd for i in eachindex(x)
        tot += x[i]
    end
    return tot
end


function _apply_fluctuations(t::VLBIImagePriors.LogRatioTransform, mimg::AbstractArray, δ::AbstractArray)
    @argcheck _checknorm(mimg) "Mean image must have unit flux when using log-ratio transformations in apply_fluctuations"
    r = to_simplex(t, baseimage(δ))
    r .= r.*baseimage(mimg)
    r .= r./_fastsum(r)
    return r
end
