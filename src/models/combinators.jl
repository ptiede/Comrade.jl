export add, convolved, smoothed, components

"""
    $(TYPEDEF)
Abstract type that denotes a composite model. Where we have
combined two models together.

# Implementation
Any implementation of a composite type must define the following methods:

- visibility_point
- uv_combinator
- imanalytic
- visanalytic
- ComradeBase.intensity_point if model intensity is `IsAnalytic`
- intensitymap! if model intensity is `NotAnalytic`
- flux

"""
abstract type CompositeModel{M1,M2} <: AbstractModel end


function modelimage(::NotAnalytic,
    model::CompositeModel,
    image; alg=FFT())

    m1 = @set model.m1 = modelimage(model.m1, image; alg)
    @set m1.m2 = modelimage(m1.m2, image; alg)
end

radialextent(m::CompositeModel) = max(radialextent(m.m1), radialextent(m.m2))

@inline visanalytic(::Type{<:CompositeModel{M1,M2}}) where {M1,M2} = visanalytic(M1)*visanalytic(M2)
@inline imanalytic(::Type{<:CompositeModel{M1,M2}}) where {M1,M2} = imanalytic(M1)*imanalytic(M2)


"""
    $(TYPEDEF)
Adds two models together to create composite models. Note that
I may change this in the future so make it easier on the compiler,
i.e. make the composite model a fancy model vector and heap allocate
stuff. This should help when combining multiple models together.
"""
struct AddModel{T1,T2} <: CompositeModel{T1,T2}
    m1::T1
    m2::T2
end
Base.:+(m1::T1, m2::T2) where {T1<:AbstractModel, T2<:AbstractModel} = AddModel(m1, m2)
Base.:-(m1, m2) = AddModel(m1, -1.0*m2)
add(m1::M1, m2::M2) where {M1<:AbstractModel, M2<:AbstractModel} = AddModel(m1, m2)

@inline combinator(::Type{<:AddModel}) = add


"""
    $(SIGNATURES)
Returns the components for a composite model. This
will return a Tuple with all the models you have constructed.
"""
components(m::AbstractModel) = (m,)
components(m::CompositeModel{M1,M2}) where
    {M1<:AbstractModel, M2<:AbstractModel} = (components(m.m1)..., components(m.m2)...)

flux(m::AddModel) = flux(m.m1) + flux(m.m2)

function intensity(m::AddModel{T1,T2}, x, y, args...) where {T1, T2}
    return intensity(m.m1, x, y, args...) + intensity(m.m2, x, y, args...)
end

function intensitymap(m::AddModel, fovx::Real, fovy::Real, nx::Int, ny::Int; pulse=DeltaPulse())
    sim1 = intensitymap(m.m1, fovx, fovy, nx, ny; pulse)
    sim2 = intensitymap(m.m2, fovx, fovy, nx, ny; pulse)
    return sim1 .+ sim2
end

function intensitymap!(sim::IntensityMap, m::AddModel)
    csim = deepcopy(sim)
    intensitymap!(csim, m.m1)
    sim .= csim
    intensitymap!(csim, m.m2)
    sim .= sim .+ csim
    return sim
end

@inline uv_combinator(::AddModel) = Base.:+

# @inline function _visibilities(model::CompositeModel{M1,M2}, u, v, t, ν, cache) where {M1,M2}
#     _combinatorvis(visanalytic(M1), visanalytic(M2), uv_combinator(model), model, u, v, t, ν, cache)
# end

@inline function visibility_point(model::CompositeModel{M1,M2}, u, v, args...) where {M1,M2}
    f = uv_combinator(model)
    v1 = visibility(model.m1, u, v, args...)
    v2 = visibility(model.m2, u, v, args...)
    return f(v1,v2)
end

@inline function intensity_point(model::CompositeModel, u, v)
    f = uv_combinator(model)
    v1 = intensity_point(model.m1, u, v)
    v2 = intensity_point(model.m2, u, v)
    return f(v1,v2)
end



"""
    $(TYPEDEF)
Smooths a `model` with a Gaussian kernel with standard deviation
of `σ`.

## Notes

```julia
m = smoothed(Disk(), 5.0)
intensity(m, 2.0, 2.0; fov=10.0, npix=128)
```
will compute the intensity of the smoothed disk with using an interpolation
with `128` pixel nodes and a total `fov` of 10 in x and y direction.

**This needs to be improved**
"""
struct ConvolvedModel{M1, M2} <: CompositeModel{M1,M2}
    m1::M1
    m2::M2
end
convolved(m1, m2) = ConvolvedModel(m1, m2)

"""
    $(SIGNATURES)
Smooths a model `m` with a Gaussian kernel with standard deviation `σ`.

# Notes
This will created a convolved model under the hood.
"""
smoothed(m, σ::Number) = convolved(m, stretched(Gaussian(), σ, σ))

@inline imanalytic(::Type{<:ConvolvedModel}) = NotAnalytic()


@inline uv_combinator(::ConvolvedModel) = Base.:*

flux(m::ConvolvedModel) = flux(m.m1)*flux(m.m2)

function intensitymap(model::ConvolvedModel, fovx::Real, fovy::Real, nx::Int, ny::Int; pulse=DeltaPulse())
    T = typeof(visibility(model, 0.0, 0.0))
    vis1 = fouriermap(model.m1, fovx, fovy, nx, ny)
    vis2 = fouriermap(model.m2, fovx, fovy, nx, ny)
    vis = ifftshift(phasedecenter!(vis1.*vis2, fovx, fovy, nx, ny))
    img = ifft(vis)
    return IntensityMap(real.(img)./(nx*ny), fovx, fovy, pulse)
end

function intensitymap!(sim::IntensityMap, model::ConvolvedModel)
    ny, nx = size(sim)
    vis1 = fouriermap(model.m1, sim)
    vis2 = fouriermap(model.m2, sim)
    vis = ifftshift(phasedecenter!(vis1.*vis2, fovx, fovy, nx, ny))
    ifft!(vis)
    for I in CartesianIndices(sim)
        sim[I] = real(vis[I])/(nx*ny)
    end
end


# function _combinatorvis(::IsAnalytic, ::IsAnalytic, f::F, m, u, v, t, ν, cache) where {F}
#     @assert length(u) == length(u) "Number of visibilities must equal number of points"
#     ff(u, v, t, ν) = f(visibility_point(m.m1, u, v, t, ν, cache), visibility_point(m.m2, u, v, t, ν, cache))
#     return mappedarray(ff, u, v, t, ν)
# end

# function _combinatorvis(::IsAnalytic, ::NotAnalytic, f::F, m, u, v, args...) where {F}
#     vis = _visibilities(m.m1, u, v, args...)
#     for i in eachindex(vis)
#         vis[i] = f(visibility_point(m.m2, u[i], v[i], args...), vis[i])
#     end
# end

# function _combinatorvis(::NotAnalytic, ::IsAnalytic, f::F, m, u, v, args...) where {F}
#     vis = _visibilities(m.m2, u, v, args...)
#     for i in eachindex(vis)
#         vis[i] = f(vis[i], visibility_point(m.m1, u[i], v[i], args...))
#     end
# end

# function _combinatorvis(::NotAnalytic, ::NotAnalytic, m, f::F, u, v, args...) where {F}
#     vis1 = _visibilities(m.m1, u, v, args...)
#     vis2 = _visibilities(m.m2, u, v, args...)
#     vis1 .= f.(vis1, vis2)
# end
