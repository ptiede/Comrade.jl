export ThreadedModel

"""
    ThreadedModel
Experimental model wrapper than enables multi-threading when evaluating `intensitymap`
"""
struct ThreadedModel{M} <: AbstractModel
    model::M
end

@inline ThreadedModel(model::ThreadedModel) = model
@inline ThreadedModel(model::M) where {M<:CompositeModel} = ThreadedModel(imanalytic(M), model)
@inline ThreadedModel(::IsAnalytic, model::CompositeModel)        = ThreadedModel{typeof(model)}(model)
@inline ThreadedModel(::NotAnalytic, model::AddModel)             = AddModel(ThreadedModel(model.m1), ThreadedModel(model.m2))
@inline ThreadedModel(::NotAnalytic, model::ConvolvedModel)       = ConvolvedModel(ThreadedModel(model.m1), ThreadedModel(model.m2))
@inline ThreadedModel(model::ModelImage) = @set model.model = ThreadedModel(model.model)

Base.@constprop :aggressive @inline visanalytic(::Type{<:ThreadedModel{M}}) where {M} = visanalytic(M)
Base.@constprop :aggressive @inline imanalytic(::Type{<:ThreadedModel{M}}) where {M} = imanalytic(M)
Base.@constprop :aggressive @inline ispolarized(::Type{<:ThreadedModel{M}}) where {M} = ispolarized(M)


@inline visibility_point(m::ThreadedModel, u, v, time, freq) = visibility_point(m.model, u, v, time, freq)
@inline intensity_point(m::ThreadedModel, p) = intensity_point(m.model, p)

@inline radialextent(m::ThreadedModel) = radialextent(basemodel(m))
@inline flux(m::ThreadedModel) = flux(basemodel(m))

using AxisKeys: keyless_unname
function intensitymap(::IsAnalytic, s::ThreadedModel, g::GriddedKeys)
    T = typeof(intensity_point(s, (X=g.X[begin], Y=g.Y[begin])))
    img = IntensityMap(Array{T}(undef, length(g.X), length(g.Y)), g)
    return intensitymap!(IsAnalytic(), img, s)
end

intensitymap(m, p, threaded::Bool) = intensitymap(m, p, static(threaded))
intensitymap(m, p, ::False) = intensitymap(m, p)
intensitymap(m, p, ::True)  = intensitymap(ThreadedModel(m), p)

intensitymap!(img::IntensityMapTypes, m, threaded::Bool) = intensitymap!(img, m, static(threaded))
intensitymap!(img::IntensityMapTypes, m, ::False) = intensitymap!(img, m)
intensitymap!(img::IntensityMapTypes, m, ::True)  = intensitymap!(img, ThreadedModel(m))

function intensitymap!(::IsAnalytic, img::IntensityMap, s::ThreadedModel)
    dx, dy = pixelsizes(img)
    mm = Base.Fix1(intensity_point, s)
    g = imagegrid(img)
    Threads.@threads for I in CartesianIndices(img)
        @inbounds img[I] = mm(g[I])*dx*dy
    end
    return img
end

function fouriermap(::IsAnalytic, m::ThreadedModel, dims::AbstractDims)
    X = dims.X
    Y = dims.Y
    uu,vv = uviterator(length(X), step(X), length(Y), step(Y))
    uvgrid = ComradeBase.grid(U=uu, V=vv)
    T = typeof(visibility(m, uvgrid[1]))
    vis = similar(uvgrid, T)
    Threads.@threads for I in CartesianIndices(vis)
        vis[I] = visibility(m, uvgrid[I])
    end
    return vis
end
