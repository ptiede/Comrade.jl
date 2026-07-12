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

# ----- flat (TransformVariables) instrument transforms ---------------------

abstract type AbstractInstrumentTransform <: TV.VectorTransform end

function TV.transform_with(flag::TV.LogJacFlag, m::AbstractInstrumentTransform, x, index)
    y, ℓ, index = _instrument_transform_with(flag, m, x, index)
    sm = site_map(m)
    sa = SiteArray(y, sm)
    return sa, ℓ, index
end

EnzymeRules.inactive_type(::Type{AbstractArray{<:IntegrationTime}}) = true
EnzymeRules.inactive_type(::Type{AbstractArray{<:FrequencyChannel}}) = true

function TV.inverse_at!(x::AbstractArray, index, t::AbstractInstrumentTransform, y::SiteArray)
    itrf = inner_transform(t)
    return TV.inverse_at!(x, index, itrf, parent(y))
end


struct InstrumentTransform{T, L <: SiteLookup} <: AbstractInstrumentTransform
    inner_transform::T
    site_map::L
end

function TV.inverse_eltype(::AbstractInstrumentTransform, x::Type)
    return eltype(x)
end


struct MarkovInstrumentTransform{T, L <: SiteLookup} <: AbstractInstrumentTransform
    inner_transform::T
    site_map::L
end

TV.dimension(m::AbstractInstrumentTransform) = TV.dimension(inner_transform(m))


@inline function _instrument_transform_with(flag::TV.LogJacFlag, m::InstrumentTransform, x, index)
    itrf = inner_transform(m)
    return TV.transform_with(flag, itrf, x, index)
end


function branchcut(x::T) where {T}
    xmod = mod2pi(x)
    return xmod > π ? xmod - convert(T, 2π) : xmod
end

@inline function _instrument_transform_with(flag::TV.LogJacFlag, m::MarkovInstrumentTransform, x, index)
    y, ℓ, index = TV.transform_with(flag, inner_transform(m), x, index)
    return markov_accumulate!(y, site_map(m)), ℓ, index
end


# ----- Std-space (ProbabilityTransports) instrument transforms -------------

abstract type AbstractStdInstrumentTransform <: PT.AbstractTransport end

site_map(t::Union{AbstractInstrumentTransform, AbstractStdInstrumentTransform}) = t.site_map
EnzymeRules.inactive(::typeof(site_map), args...) = nothing
inner_transform(t::Union{AbstractInstrumentTransform, AbstractStdInstrumentTransform}) = t.inner_transform

PT.dimension(m::AbstractStdInstrumentTransform) = PT.dimension(inner_transform(m))
PT.pback_eltype(m::AbstractStdInstrumentTransform) =
    PT.pback_eltype(inner_transform(m))


struct StdInstrumentTransform{T, L <: SiteLookup} <: AbstractStdInstrumentTransform
    inner_transform::T
    site_map::L
end

struct StdMarkovInstrumentTransform{T, L <: SiteLookup} <: AbstractStdInstrumentTransform
    inner_transform::T
    site_map::L
end

function PT.pfwd_step(m::StdInstrumentTransform, x, index)
    y, index = PT.pfwd_step(inner_transform(m), x, index)
    return SiteArray(y, site_map(m)), index
end

function PT.pfwd_step(m::StdMarkovInstrumentTransform, x, index)
    y, index = PT.pfwd_step(inner_transform(m), x, index)
    return SiteArray(markov_accumulate!(y, site_map(m)), site_map(m)), index
end

function PT.pback_step!(y::AbstractVector, index, m::StdInstrumentTransform, x::SiteArray)
    return PT.pback_step!(y, index, inner_transform(m), parent(x))
end

function PT.pback_step!(y::AbstractVector, index, m::StdMarkovInstrumentTransform, x::SiteArray)
    yd = markov_difference(x, site_map(m))
    return PT.pback_step!(y, index, inner_transform(m), parent(yd))
end


# ----- shared site-mapping helpers ----------------------------------------

EnzymeRules.inactive_type(::Type{<:SiteLookup}) = true

# Cumulative-sum each site's parameters (the Markov structure) and wrap phases to (-π, π].
@inline function markov_accumulate!(y, site_map::SiteLookup)
    yout = site_sum(y, site_map)
    yout .= branchcut.(yout)
    return yout
end

# The inverse of `markov_accumulate!`'s cumulative sum: recover the raw per-site increments.
# One `copy` guards the caller's input; `site_diff!` then differences in place (no per-site
# temporaries).
function markov_difference(y, site_map::SiteLookup)
    yd = copy(y)
    site_diff!(yd, site_map)
    return yd
end

@inline function site_sum(y, site_map::SiteLookup)
    vals = values(lookup(site_map))
    @inbounds for site in vals
        acc = zero(eltype(y))
        @inbounds ys = @view y[site]
        for i in eachindex(ys)
            acc += ys[i]
            ys[i] = acc
        end
    end
    return y
end


function site_diff!(y, site_map::SiteLookup)
    map(site_map.lookup) do site
        simplediff!(@view y[site])
    end
    return nothing
end

# In-place first difference: `x[i] -> x[i] - x[i-1]`, leaving `x[begin]` untouched. Iterating
# from the end keeps each `x[i-1]` at its original value when it is used, so no copy is needed.
function simplediff!(x::AbstractVector)
    @inbounds for i in reverse((firstindex(x) + 1):lastindex(x))
        x[i] -= x[i - 1]
    end
    return x
end

function TV.inverse_at!(x::AbstractArray, index, t::MarkovInstrumentTransform, y::SiteArray)
    yd = markov_difference(y, site_map(t))
    return TV.inverse_at!(x, index, inner_transform(t), yd)
end
