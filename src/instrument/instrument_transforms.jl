# Instrument-model parameter transforms.
#
# These map between a flat parameter vector and a `SiteArray` of instrument
# parameters. There are two families, one per latent space:
#
#   * the flat (`TVFlat`/NUTS) path is a `TV.VectorTransform` so the whole flat
#     prior tree is a single TransformVariables transform (which threads the
#     Jacobian); and
#   * the Std-space (`StdUniform`/`StdNormal`, e.g. nested sampling) path is a
#     `ProbabilityTransports.AbstractTransport` node, whose `pfwd_step` /
#     `pback_step!` mirror the flat `transform_with` / `inverse_at!` but carry
#     no Jacobian (the transport is exact in those spaces).

# ----- flat (TransformVariables) instrument transform ----------------------

struct InstrumentTransform{T, L <: SiteLookup} <: TV.VectorTransform
    inner_transform::T
    site_map::L
end

EnzymeRules.inactive_type(::Type{AbstractArray{<:IntegrationTime}}) = true
EnzymeRules.inactive_type(::Type{AbstractArray{<:FrequencyChannel}}) = true

TV.dimension(m::InstrumentTransform) = TV.dimension(inner_transform(m))
TV.inverse_eltype(::InstrumentTransform, x::Type) = eltype(x)

function TV.transform_with(flag::TV.LogJacFlag, m::InstrumentTransform, x, index)
    y, ℓ, index = TV.transform_with(flag, inner_transform(m), x, index)
    return SiteArray(y, site_map(m)), ℓ, index
end

function TV.inverse_at!(x::AbstractArray, index, t::InstrumentTransform, y::SiteArray)
    return TV.inverse_at!(x, index, inner_transform(t), parent(y))
end


# ----- Std-space (ProbabilityTransports) instrument transform ---------------

struct StdInstrumentTransform{T, L <: SiteLookup} <: PT.AbstractTransport
    inner_transform::T
    site_map::L
end

PT.dimension(m::StdInstrumentTransform) = PT.dimension(inner_transform(m))
PT.pback_eltype(m::StdInstrumentTransform) = PT.pback_eltype(inner_transform(m))

function PT.pfwd_step(m::StdInstrumentTransform, x, index)
    y, index = PT.pfwd_step(inner_transform(m), x, index)
    return SiteArray(y, site_map(m)), index
end

function PT.pback_step!(y::AbstractVector, index, m::StdInstrumentTransform, x::SiteArray)
    return PT.pback_step!(y, index, inner_transform(m), parent(x))
end


site_map(t::Union{InstrumentTransform, StdInstrumentTransform}) = t.site_map
EnzymeRules.inactive(::typeof(site_map), args...) = nothing
inner_transform(t::Union{InstrumentTransform, StdInstrumentTransform}) = t.inner_transform

EnzymeRules.inactive_type(::Type{<:SiteLookup}) = true
