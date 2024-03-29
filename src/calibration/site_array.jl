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
    return SiteArray(getindex(parent(m), I...), getindex(m.time, I...), getindex(m.sites, I...))
end

function Base.view(A::SiteArray, I...)
    return SiteArray(view(A.data, I...), view(A.time, I...), view(A.sites, I...))
end


function Base.similar(m::SiteArray, ::Type{S}) where {S}
    return SiteArray(similar(parent(m), S), m.time, m.sites)
end

function Base.similar(m::SiteArray, ::Type{S}, dims::Dims) where {S}
    any(x->x[1]!=x[2], zip(dims, size(m))) && throw(DimensionMismatch("Size of new array must be a identical to passed SiteArray"))
    return SiteArray(similar(parent(m), S, dims), m.time, m.sites)
end



Base.BroadcastStyle(::Type{<:SiteArray}) = Broadcast.ArrayStyle{SiteArray}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SiteArray}}, ::Type{ElType}) where {ElType}
    # Scan inputs for the time and sites
    sarr = find_sitesarr(bc)
    return SiteArray(similar(parent(sarr), ElType), sarr.time, sarr.sites)
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


function time_interval(arr::SiteArray, a::Interval)
    inds = findall(in(a), times(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end


function Base.getindex(arr::SiteArray; site=Base.Colon(), time=Base.Colon())
    if site == Base.Colon() && time == Base.Colon()
        return arr
    end
    if site == Base.Colon()
        return time_interval(arr, time)
    end
    if time == Base.Colon()
        return Comrade.site(arr, site)
    end

    return site_time_interval(arr, site, time)
end

function site_time_interval(arr::SiteArray, site::Union{Symbol, String}, time::Interval)
    site_time_interval(arr, (site,), time)
end


function site_time_interval(arr::SiteArray, site, time::Interval)
    inds = findall(i->((Comrade.sites(arr)[i] ∈ site)&&(Comrade.times(arr)[i] ∈ time)), eachindex(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end

struct SiteMap{L<:NamedTuple, N, Ti<:AbstractArray{<:Number, N}, Fr<:AbstractArray{<:Number, N}, Sy<:AbstractArray{<:Any, N}}
    lookup::L
    times::Ti
    frequencies::Fr
    sites::Sy
end

function sitemap!(f, out::AbstractArray, gains::AbstractArray, cache::AbstractJonesCache)
    sm = cache.site_map
    map(sm.lookup) do site
        ys = @view gains[site]
        outs = @view out[site]
        outs .= f.(ys)
    end
end

function sitemap(f, gains::AbstractArray, cache::AbstractJonesCache)
    out = similar(T, gains)
    sitemap!(f, out, gains, cache)
    return out
end

function sitemap!(::typeof(cumsum), out::AbstractArray, gains::AbstractArray, cache::AbstractJonesCache)
    map(site_map.lookup) do site
        ys = @view gains[site]
        cumsum!(ys, ys)
    end
    return out
end

function sitemap(s::SiteArray)
    slist = Comrade.sites(s)
    sites = Tuple(sort(unique(slist)))
    return SiteMap(NamedTuple{sites}(map(p->findall(==(p), slist), sites)), times(s), frequencies(s), slist)
end
