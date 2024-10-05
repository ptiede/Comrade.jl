abstract type AbstractInstrumentTransform <: TV.VectorTransform end
site_map(t::AbstractInstrumentTransform) = t.site_map
inner_transform(t::AbstractInstrumentTransform) = t.inner_transform

function TV.transform_with(flag::TV.LogJacFlag, m::AbstractInstrumentTransform, x, index)
    y, ℓ, index = _instrument_transform_with(flag, m, x, index)
    sm = m.site_map
    return SiteArray(y, sm.times, sm.frequencies, sm.sites), ℓ, index
end

function TV.inverse_at!(x::AbstractArray, index, t::AbstractInstrumentTransform, y::SiteArray)
    itrf = inner_transform(t)
    return TV.inverse_at!(x, index, itrf, parent(y))
end


struct InstrumentTransform{T, L<:SiteLookup} <: AbstractInstrumentTransform
    inner_transform::T
    site_map::L
end

function TV.inverse_eltype(::AbstractInstrumentTransform, x)
    return eltype(x)
end


struct MarkovInstrumentTransform{T, L<:SiteLookup} <: AbstractInstrumentTransform
    inner_transform::T
    site_map::L
end

TV.dimension(m::AbstractInstrumentTransform) = TV.dimension(inner_transform(m))


@inline function _instrument_transform_with(flag::TV.LogJacFlag, m::InstrumentTransform, x, index)
    itrf = inner_transform(m)
    return TV.transform_with(flag, itrf, x, index)
end

@inline function _instrument_transform_with(flag::TV.LogJacFlag, m::MarkovInstrumentTransform, x, index)
    (;inner_transform, site_map) = m
    y, ℓ, index = TV.transform_with(flag, inner_transform, x, index)
    yout = site_sum(y, site_map)
    return yout, ℓ, index
end

@inline function site_sum(y, site_map::SiteLookup)
    yout = similar(y)
    @inbounds for site in site_map.lookup
        ys = @view y[site]
        # y should never alias so we should be fine here.
        youts = @view yout[site]
        cumsum!(youts, ys)
    end
    return yout
end

# function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(_instrument_transform_with), flag, m::MarkovInstrumentTransform, x, index)
#     (y, ℓ, index2), dt = rrule_via_ad(config, TV.transform_with, flag, m.inner_transform, x, index)
#     site_sum!(y, m.site_map)
#     py = ProjectTo(y)
#     function _markov_transform_pullback(Δ)
#         Δy  = similar(y)
#         Δy .= py(unthunk(Δ[1]))
#         autodiff(Reverse, site_sum!, Const, Duplicated(y, Δy), Const(m.site_map))
#         din = dt((Δy, Δ[2], NoTangent()))
#         return din
#     end
#     return (y, ℓ, index2), _markov_transform_pullback
# end


function site_diff!(y, site_map::SiteLookup)
    map(site_map.lookup) do site
        ys = @view y[site]
        simplediff!(ys, copy(ys))
    end
    return nothing
end

function simplediff(x::AbstractVector)
    y = zero(x)
    simplediff!(y, x)
    return y
end

function simplediff!(y::AbstractVector, x::AbstractVector)
    y[begin] = x[begin]
    y[begin+1:end] .= @views x[begin+1:end] .- @views x[begin:end-1]
    return nothing
end

function TV.inverse_at!(x::AbstractArray, index, t::MarkovInstrumentTransform, y::SiteArray)
    (;inner_transform, site_map) = t
    # Now difference to get the raw values
    yd = copy(y)
    site_diff!(yd, site_map)
    # and now inverse
    return TV.inverse_at!(x, index, inner_transform, yd)
end
