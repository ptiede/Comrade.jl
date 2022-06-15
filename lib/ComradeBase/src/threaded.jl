export ThreadedModel

"""
    ThreadedModel
Experimental model wrapper than enables multi-threading. This will likely
be depreciated in the future for moving threading to executors similar to FLoops.jl
"""
struct ThreadedModel{M} <: AbstractModel
    model::M
end
@inline imanalytic(::Type{ThreadedModel{M}}) where {M} = imanalytic(M)
@inline visanalytic(::Type{ThreadedModel{M}}) where {M} = visanalytic(M)
@inline isprimitive(::Type{ThreadedModel{M}}) where {M} = isprimitive(M)


@inline visibility_point(m::ThreadedModel, u, v) = visibility_point(m.model, u, v)
@inline intensity_point(m::ThreadedModel, x, y) = intensity_point(m.model, x, y)


function intensitymap(::IsAnalytic, s::ThreadedModel, fovx::Number, fovy::Number, nx::Int, ny::Int; executor= pulse=ComradeBase.DeltaPulse())
    T = typeof(intensity_point(s, 0.0, 0.0))
    img = IntensityMap(zeros(T, ny, nx), fovx, fovy, pulse)
    intensitymap!(IsAnalytic(), img, s)
    return img
end

function intensitymap!(::IsAnalytic, img::AbstractIntensityMap, s::ThreadedModel)
    x,y = imagepixels(img)
    dx, dy = pixelsizes(img)
    Threads.@threads for I in CartesianIndices(img)
        iy,ix = Tuple(I)
        img[I] = intensity_point(s, x[ix], y[iy])*dx*dy
    end
    return img
end
