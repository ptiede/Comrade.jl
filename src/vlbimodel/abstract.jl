abstract type AbstractVLBIModel{M, J} <: AbstractModel end
skymodel(m::AbstractVLBIModel) = getfield(m, :sky)
skymodel(m::AbstractVLBIModel, p) = skymodel(m)(p)
flux(m::AbstractVLBIModel) = flux(basemodel(m))
radialextent(m::AbstractVLBIModel) = radialextent(basemodel(m))
instrumentmodel(m::AbstractVLBIModel) = getfield(m, :instrument)
instrumentmodel(m::AbstractVLBIModel, p) = instrumentmodel(m)(p)

visanalytic(::Type{<:AbstractVLBIModel{J,M}}) where {J,M} = visanalytic(M)
imanalytic(::Type{<:AbstractVLBIModel{J,M}}) where {J,M} = imanalytic(M)
ispolarized(::Type{<:AbstractVLBIModel{J,M}}) where {J,M} = ispolarized(M)


struct VLBIModel{S<:SkyModel, J<:InstrumentModel, D}
    sky::S
    instrument::J
end

function forward_model(m::VLBIModel, params)
    skym = skymodel(m, params.sky)
    vis = visibilitymap(skym, domain)
    return apply_instrument(vis, instrumentmodel(m), params.instrument)
end


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
