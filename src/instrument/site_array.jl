export SiteArray, site, time_interval

"""
    SiteArray(data, times, frequencies, sites)

A `SiteArray` is an array of data that has a specified ordering of `times`, `frequencies`, and `sites`.
Each data point is assigned a unique `time`, `frequency`, and `site` code. This allows for easy
selection of data points based on these criteria and forms the base array for instrument modeling.

To select a subset of the data based on a specifid site, time and frequency you can use
```
sarr[S=:ALMA, Ti=1:10, Fr=1:10]
```
which will grab the first 10 time and frequency points for the ALMA site.

Otherwise indexing into the array will return an element whose time, frequency, and site are
the element of the `times`, `frequencies`, and `sites` arrays respectively.
"""
struct SiteArray{T, N, A<:AbstractArray{T,N}, Ti<:AbstractArray{<:IntegrationTime, N}, Fr<:AbstractArray{<:FrequencyChannel, N}, Sy<:AbstractArray{<:Any, N}} <: AbstractArray{T, N}
    data::A
    times::Ti
    frequencies::Fr
    sites::Sy
end

# function ChainRulesCore.rrule(::Type{SiteArray}, data::AbstractArray, args...)
#     s = SiteArray(data, args...)
#     pd = ProjectTo(data)
#     function _SiteArrayPB(Δ)
#         # @info typeof(Δ)
#         (NoTangent(), @thunk(pd(Δ)), map(i->NoTangent(), args)...)
#     end
#     return s, _SiteArrayPB
# end

times(a::SiteArray) = a.times
sites(a::SiteArray) = a.sites
frequencies(a::SiteArray) = a.frequencies

EnzymeRules.inactive(::(typeof(Base.size)), ::SiteArray) = nothing
Base.parent(a::SiteArray) = getfield(a, :data)
Base.size(a::SiteArray) = size(parent(a))
Base.IndexStyle(::Type{<:SiteArray{T, N, A}}) where {T, N, A} = Base.IndexStyle(A)
Base.@propagate_inbounds Base.getindex(a::SiteArray{T}, i::Integer) where {T} = getindex(parent(a), i)
Base.@propagate_inbounds Base.getindex(a::SiteArray, I::Vararg{Integer, N}) where {N} = getindex(parent(a), I...)
Base.@propagate_inbounds Base.setindex!(m::SiteArray, v, i::Integer) = setindex!(parent(m), v, i)
Base.@propagate_inbounds Base.setindex!(m::SiteArray, v, i::Vararg{Integer, N}) where {N} = setindex!(parent(m), v, i...)
Base.@propagate_inbounds function Base.getindex(m::SiteArray, I...)
    return SiteArray(getindex(parent(m), I...), getindex(m.times, I...), getindex(m.frequencies, I...), getindex(m.sites, I...))
end

Base.@propagate_inbounds function Base.view(A::SiteArray, I...)
    return SiteArray(view(A.data, I...), view(times(A), I...), view(frequencies(A), I...), view(sites(A), I...))
end

# function ChainRulesCore.ProjectTo(s::SiteArray)
#     return ProjectTo{SiteArray}(; data=parent(s),
#                                   times=times(s),
#                                   frequencies=frequencies(s),
#                                   sites=sites(s))
# end

# (project::ProjectTo{SiteArray})(s) = SiteArray(s, project.times, project.frequencies, project.sites)
# (project::ProjectTo{SiteArray})(s::SiteArray) = s
# (project::ProjectTo{SiteArray})(s::AbstractZero) = s
# (project::ProjectTo{SiteArray})(s::Tangent) = SiteArray(s.data, project.times, project.frequencies, project.sites)


EnzymeRules.inactive(::typeof(times), ::SiteArray) = nothing
EnzymeRules.inactive(::typeof(frequencies), ::SiteArray) = nothing
EnzymeRules.inactive(::typeof(sites), ::SiteArray) = nothing

# ntzero(x::NamedTuple) = map(ntzero, x)
# ntzero(x::Tuple) = map(ntzero, x)
# ntzero(x) = zero(x)


function Base.similar(m::SiteArray, ::Type{S}) where {S}
    return SiteArray(similar(parent(m), S), m.times, m.frequencies, m.sites)
end

function Base.similar(m::SiteArray, ::Type{S}, dims::Dims) where {S}
    any(x->x[1]!=x[2], zip(dims, size(m))) && throw(DimensionMismatch("Size of new array must be a identical to passed SiteArray"))
    return SiteArray(similar(parent(m), S, dims), m.times, m.frequencies, m.sites)
end



Base.BroadcastStyle(::Type{<:SiteArray}) = Broadcast.ArrayStyle{SiteArray}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SiteArray}}, ::Type{ElType}) where {ElType}
    # Scan inputs for the time and sites
    sarr = find_sitesarr(bc)
    return SiteArray(similar(parent(sarr), ElType), sarr.times, sarr.frequencies, sarr.sites)
end

find_sitesarr(bc::Broadcast.Broadcasted) = find_sitesarr(bc.args)
find_sitesarr(args::Tuple) = find_sitesarr(find_sitesarr(args[1]), Base.tail(args))
find_sitesarr(x) = x
find_sitesarr(::Tuple{}) = nothing
find_sitesarr(x::SiteArray, rest) = x
find_sitesarr(::Any, rest) = find_sitesarr(rest)

function site(arr::SiteArray, p::Union{Symbol, String})
    return site(arr, (p,))
end

function site(arr::SiteArray, p)
    inds = findall(in(p), sites(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end


function time(arr::SiteArray, a::Union{AbstractInterval, Real})
    inds = findall(in(a), times(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end

function frequency(arr::SiteArray, a::Union{AbstractInterval, Real})
    inds = findall(in(a), times(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end

_equalorin(x::T, y::T) where {T} = x == y 
_equalorin(x::Real, y) = x ∈ y
_equalorin(x, y::Real) = y ∈ x
_equalorin(x, y) = y ∈ x
_equalorin(x, ::typeof(Base.Colon())) = true
_equalorin(::typeof(Base.Colon()), x) = true
const Indexable = Union{Integer, AbstractArray{<:Integer}, BitArray}

function Base.getindex(arr::SiteArray; Fr=Base.Colon(), S=Base.Colon(), Ti=Base.Colon())
    Fr2 = isa(Fr, Indexable) ? unique(arr.frequencies)[Fr] : Fr
    S2 = isa(S, Indexable) ? unique(arr.sites)[S] : S
    Ti2 = isa(Ti, Indexable) ? unique(arr.times)[Ti] : Ti
    inds = findall(i->(_equalorin(S2, Comrade.sites(arr)[i])&&_equalorin(Ti2, Comrade.times(arr)[i])&&_equalorin(Fr2, Comrade.frequencies(arr)[i])), eachindex(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end


struct SiteLookup{L<:NamedTuple, N, Ti<:AbstractArray{<:IntegrationTime, N}, Fr<:AbstractArray{<:FrequencyChannel, N}, Sy<:AbstractArray{<:Any, N}}
    lookup::L
    times::Ti
    frequencies::Fr
    sites::Sy
end

function sitemap!(f, out::AbstractArray, gains::AbstractArray, slook::SiteLookup)
    map(slook.lookup) do site
        ysite = @view gains[site]
        outsite = @view out[site]
        outsite .= f.(ysite)
    end
end

function sitemap(f, gains::AbstractArray{T}, slook::SiteLookup) where {T}
    out = similar(gains)
    sitemap!(f, out, gains, slook)
    return out
end

function sitemap!(::typeof(cumsum), out::AbstractArray, gains::AbstractArray, slook::SiteLookup)
    map(slook.lookup) do site
        ys = @view gains[site]
        cumsum!(ys, ys)
    end
    return out
end

"""
    SiteLookup(s::SiteArray)

Construct a site lookup dictionary for a site array.
"""
function SiteLookup(s::SiteArray)
    return SiteLookup(times(s), frequencies(s), sites(s))
end

function SiteLookup(times::AbstractVector, frequencies::AbstractArray, sites::AbstractArray)
    slist = Tuple(sort(unique(sites)))
    return SiteLookup(NamedTuple{slist}(map(p->findall(==(p), sites), slist)), times, frequencies, sites)
end

"""
    SiteArray(arr, sitelookup::SiteLookup)

Construct a site array with the entries `arr` and the site ordering implied by
`sitelookup`.
"""
function SiteArray(a::AbstractArray, map::SiteLookup)
    return SiteArray(a, map.times, map.frequencies, map.sites)
end

function SiteArray(data::SiteArray{T, N},
                   times::AbstractArray{<:IntegrationTime, N}, 
                   frequencies::AbstractArray{<:FrequencyChannel, N}, 
                   sites::AbstractArray{<:Number, N}) where {T, N}
    return data
end

function SiteArray(data::SiteArray, ::SiteLookup)
    return data
end
