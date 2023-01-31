export ThreadedModel

"""
    ThreadedModel
Experimental model wrapper than enables multi-threading when evaluating `intensitymap`
"""
struct ThreadedModel{M} <: AbstractModel
    model::M
end


Base.@constprop :aggressive @inline visanalytic(::Type{<:AbstractModifier{M}}) where {M} = visanalytic(M)
Base.@constprop :aggressive @inline imanalytic(::Type{<:AbstractModifier{M}}) where {M} = imanalytic(M)
Base.@constprop :aggressive @inline ispolarized(::Type{<:AbstractModifier{M}}) where {M} = ispolarized(M)


@inline visibility_point(m::ThreadedModel, u, v, time, freq) = visibility_point(m.model, u, v, time, freq)
@inline intensity_point(m::ThreadedModel, p) = intensity_point(m.model, p)

@inline radialextent(m::ThreadedModel) = radialextent(basemodel(m))
@inline flux(m::ThreadedModel) = flux(basemodel(m))

function intensitymap(::IsAnalytic, s::ThreadedModel, p::GriddedKeys)
    dx = step(dims.X)
    dy = step(dims.Y)
    img = @.. true intensity_point.(Ref(s), imagegrid(p)).*dx.*dy
    return IntensityMap(AxisKeys.keyless_unname(img), dims)
    return img
end

function intensitymap!(::IsAnalytic, img::AbstractIntensityMap, s::ThreadedModel)
    dx, dy = pixelsizes(img)
    g = imagegrid(img)
    img .= intensity_point.(Ref(s), g).*dx.*dy
    return img
end

function fouriermap(m::ThreadedModel, dims::AbstractDims)
    X = dims.X
    Y = dims.Y
    uu,vv = uviterator(length(X), step(X), length(Y), step(Y))
    uvgrid = ComradeBase.grid(U=uu, V=vv)
    vis = @.. true visibility.(Ref(m), uvgrid)
end
