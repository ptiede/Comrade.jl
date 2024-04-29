abstract type AbstractVLBIModel{M, J} <: AbstractModel end
skymodel(m::AbstractVLBIModel) = getfield(m, :sky)
flux(m::AbstractVLBIModel) = flux(basemodel(m))
radialextent(m::AbstractVLBIModel) = radialextent(basemodel(m))
instrumentmodel(m::AbstractVLBIModel) = getfield(m, :instrument)

visanalytic(::Type{<:AbstractVLBIModel{J,M}}) where {J,M} = visanalytic(M)
imanalytic(::Type{<:AbstractVLBIModel{J,M}}) where {J,M} = imanalytic(M)
ispolarized(::Type{<:AbstractVLBIModel{J,M}}) where {J,M} = ispolarized(M)


function intensitymap(model::AbstractVLBIModel, dims::AbstractDomain)
    return intensitymap(skymodel(model), dims)
end

function intensity_point(model::AbstractVLBIModel, p)
    return intensity_point(skymodel(model), p)
end

function visibilitymap_analytic(model::AbstractVLBIModel, p::AbstractDomain)
    skym = model.sky(model.skyparams)
    vis = visibilitymap_analytic(skym, p::AbstractDomain)
    return apply_instrument(vis, model.instrument, model.instrumentparams)
end

function visibilitymap_numeric(model::AbstractVLBIModel, p::AbstractDomain)
    skym = model.sky
    vis = visibilitymap_numeric(skym, p::AbstractDomain)
    return apply_instrument(vis, model.instrument)
end
