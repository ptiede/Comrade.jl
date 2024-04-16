export SiteArray, site, time_interval

struct SiteArray{T, N, A<:AbstractArray{T,N}, Ti<:AbstractArray{<:Number, N}, Fr<:AbstractArray{<:Number, N}, Sy<:AbstractArray{<:Any, N}} <: AbstractArray{T, N}
    data::A
    times::Ti
    frequencies::Fr
    sites::Sy
end

times(a::SiteArray) = a.times
sites(a::SiteArray) = a.sites
frequencies(a::SiteArray) = a.frequencies


Base.parent(a::SiteArray) = a.data
Base.size(a::SiteArray) = size(parent(a))
Base.IndexStyle(::Type{<:SiteArray{T, N, A}}) where {T, N, A} = Base.IndexStyle(A)
Base.getindex(a::SiteArray, i::Integer) = getindex(parent(a), i)
Base.getindex(a::SiteArray, I::Vararg{Int, N}) where {N} = getindex(parent(a), I...)
Base.setindex!(m::SiteArray, v, i::Int) = setindex!(parent(m), v, i)
Base.setindex!(m::SiteArray, v, i::Vararg{Int, N}) where {N} = setindex!(parent(m), v, i...)
function Base.getindex(m::SiteArray, I...)
    return SiteArray(getindex(parent(m), I...), getindex(m.times, I...), getindex(m.frequencies, I...), getindex(m.sites, I...))
end

function Base.view(A::SiteArray, I...)
    return SiteArray(view(A.data, I...), view(A.times, I...), view(A.frequencies, I...), view(A.sites, I...))
end


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


function time(arr::SiteArray, a::Union{Interval, IntegrationTime})
    inds = findall(in(a), times(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end

function frequency(arr::SiteArray, a::Union{Interval, FrequencyChannel})
    inds = findall(in(a), times(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end


@inline function _maybe_all(arr, X)
    if X isa Base.Colon
        ext = extrema(arr)
        return Interval(ext[1], ext[2])
    else
        return X
    end
end

function Base.getindex(arr::SiteArray; F=Base.Colon(), S=Base.Colon(), T=Base.Colon())
    T2 = _maybe_all(times(arr), T)
    F2 = _maybe_all(frequencies(arr), F)
    S2 = S isa Base.Colon ? unique(S) : S
    return select_region(arr, S2, T2, F2)
end

function select_region(arr::SiteArray, S::Symbol, T::Union{IntegrationTime, Interval}, F::Union{FrequencyChannel, Interval})
    select_region(arr, (S,), T, F)
end


function select_region(arr::SiteArray, site, time::Interval, T::Union{IntegrationTime, Interval}, F::Union{FrequencyChannel, Interval})
    inds = findall(i->((Comrade.sites(arr)[i] ∈ site)&&(Comrade.times(arr)[i] ∈ time)), eachindex(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end

struct SiteLookup{L<:NamedTuple, N, Ti<:AbstractArray{<:Number, N}, Fr<:AbstractArray{<:Number, N}, Sy<:AbstractArray{<:Any, N}}
    lookup::L
    times::Ti
    frequencies::Fr
    sites::Sy
end

function sitemap!(f, out::AbstractArray, gains::AbstractArray, cache::SiteLookup)
    sm = cache.site_map
    map(sm.lookup) do site
        ysite = @view gains[site]
        outsite = @view out[site]
        outsite .= f.(ysite)
    end
end

function sitemap(f, gains::AbstractArray, cache::SiteLookup)
    out = similar(T, gains)
    sitemap!(f, out, gains, cache)
    return out
end

function sitemap!(::typeof(cumsum), out::AbstractArray, gains::AbstractArray, cache::SiteLookup)
    map(site_map.lookup) do site
        ys = @view gains[site]
        cumsum!(ys, ys)
    end
    return out
end

function SiteLookup(s::SiteArray)
    return SiteLookup(times(s), frequencies(s), sites(s))
end

function SiteLookup(sites::AbstractArray, times::AbstractVector,  frequencies::AbstractArray)
    slist = Tuple(sort(unique(sites)))
    println(slist)
    return SiteLookup(NamedTuple{slist}(map(p->findall(==(p), sites), slist)), times, frequencies, sites)
end

function SiteArray(a::AbstractArray, map::SiteLookup)
    return SiteArray(a, map.times, map.frequencies, map.sites)
end
