export SiteArray, site

"""
    CodedAxis(axis, code)

A lazy dimension-coded vector: element `i` is `axis[code[i]]`. The instrument metadata
(times, frequencies, sites) stores the handful of unique axis values once and a compact
integer code per point, instead of repeating full values; views/subsets slice the codes.
This is a dedicated type rather than `view(axis, code)` so the whole family can be marked
`EnzymeRules.inactive_type` (metadata is never differentiated).
"""
struct CodedAxis{T, X <: AbstractVector{T}, C <: AbstractArray{<:Integer}} <: AbstractVector{T}
    axis::X
    code::C
end

Base.size(a::CodedAxis) = size(a.code)
Base.IndexStyle(::Type{<:CodedAxis}) = Base.IndexLinear()
Base.@propagate_inbounds Base.getindex(a::CodedAxis, i::Int) = a.axis[a.code[i]]
Base.@propagate_inbounds Base.getindex(a::CodedAxis, I::AbstractArray) = CodedAxis(a.axis, a.code[I])
Base.@propagate_inbounds Base.view(a::CodedAxis, I...) = CodedAxis(a.axis, view(a.code, I...))
Base.similar(a::CodedAxis, ::Type{S}, dims::Dims) where {S} = similar(a.axis, S, dims)

EnzymeRules.inactive_type(::Type{<:CodedAxis}) = true

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
struct SiteArray{T, N, A <: AbstractArray{T, N}, Ti <: AbstractArray{<:Number, N}, Fr <: AbstractArray{<:Number, N}, Sy <: AbstractArray{<:Any, N}} <: AbstractArray{T, N}
    data::A
    times::Ti
    frequencies::Fr
    sites::Sy
end

function Adapt.parent_type(::Type{<:Comrade.SiteArray{T, N, A}}) where {T, N, A}
    return A
end

function Adapt.adapt(to, x::SiteArray)
    data = Adapt.adapt(to, parent(x))
    return SiteArray(data, x.times, x.frequencies, x.sites)
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
Base.@propagate_inbounds Base.getindex(a::SiteArray, i::Int) = rgetindex(parent(a), i)
Base.@propagate_inbounds Base.getindex(a::SiteArray{T, N}, I::Vararg{Int, N}) where {T, N} = rgetindex(parent(a), I...)
Base.@propagate_inbounds Base.setindex!(m::SiteArray, v, i::Integer) = setindex!(parent(m), v, i)
Base.@propagate_inbounds Base.setindex!(m::SiteArray, v, i::Vararg{Integer, N}) where {N} = setindex!(parent(m), v, i...)
Base.@propagate_inbounds function Base.getindex(m::SiteArray, I::AbstractArray...)
    return SiteArray(rgetindex(parent(m), I...), rgetindex(m.times, I...), rgetindex(m.frequencies, I...), rgetindex(m.sites, I...))
end

Base.@propagate_inbounds function Base.view(A::SiteArray, I::AbstractArray...)
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
    any(x -> x[1] != x[2], zip(dims, size(m))) && throw(DimensionMismatch("Size of new array must be a identical to passed SiteArray"))
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

# Does the stored metadata value match the user's selector? Directional (selector first,
# stored value second): stored values are always scalars (Real centers or site Symbols),
# selectors may be a scalar, a collection/interval, or Colon.
_selmatch(::typeof(Base.Colon()), stored) = true
_selmatch(sel::Real, stored::Real) = sel == stored
_selmatch(sel::Symbol, stored::Symbol) = sel === stored
_selmatch(sel, stored) = stored ∈ sel
const Indexable = Union{Integer, AbstractArray{<:Integer}, BitArray}

function Base.getindex(arr::SiteArray; Fr = Base.Colon(), S = Base.Colon(), Ti = Base.Colon())
    Fr2 = isa(Fr, Indexable) ? unique(arr.frequencies)[Fr] : Fr
    S2 = isa(S, Indexable) ? unique(arr.sites)[S] : S
    Ti2 = isa(Ti, Indexable) ? unique(arr.times)[Ti] : Ti
    inds = findall(i -> (_selmatch(S2, Comrade.sites(arr)[i])&&_selmatch(Ti2, Comrade.times(arr)[i])&&_selmatch(Fr2, Comrade.frequencies(arr)[i])), eachindex(arr))
    nd = view(parent(arr), inds)
    return SiteArray(nd, view(times(arr), inds), view(frequencies(arr), inds), view(sites(arr), inds))
end


"""
    SiteFreqChain

The observed points of one (site, frequency-channel) pair: the flat indices into the
parameter storage (`inds`, ascending in time), the codes of each point into the parent
lookup's time axis (`tcodes`, ascending), and the code of the chain's frequency channel
(`fcode`). Only observed (time, frequency) combinations are represented — sites that
join/leave the array or observe a subset of the channels contribute exactly the points
they observed. Use [`chaintimes`](@ref Comrade.chaintimes)/[`chainfreq`](@ref
Comrade.chainfreq) to read the axis values.
"""
struct SiteFreqChain{V <: AbstractVector{<:Integer}, C <: AbstractVector{<:Integer}, I <: Integer}
    inds::V
    tcodes::C
    fcode::I
end

Base.length(c::SiteFreqChain) = length(c.inds)

"""
    SiteLookup

The instrument-parameter metadata in dimension-coded form: the unique time/frequency/site
axis values (with the segment edges of the time and frequency axes, used only when
mapping data onto parameter slots), compact per-point integer codes into those axes, and
the per-(site, frequency-channel) chain structure ([`sitechains`](@ref
Comrade.sitechains)). Time centers are in the same units as `array[:Ti]` (UTC hours) and
frequency centers in Hz; no per-point segment objects are stored.
"""
struct SiteLookup{L <: NamedTuple, N, TA <: AbstractVector{<:Real}, FA <: AbstractVector{<:Real}, SA <: AbstractVector{Symbol}, C <: AbstractArray{<:Integer, N}}
    chains::L
    taxis::TA  # unique time centers (sorted)
    tlo::TA    # segment edges per time-axis entry (construction/search only)
    thi::TA
    faxis::FA  # unique frequency centers (sorted)
    flo::FA
    fhi::FA
    saxis::SA  # unique sites (sorted)
    tcode::C   # per-point codes into the axes
    fcode::C
    scode::C
end

times(s::SiteLookup) = CodedAxis(s.taxis, s.tcode)
sites(s::SiteLookup) = CodedAxis(s.saxis, s.scode)
frequencies(s::SiteLookup) = CodedAxis(s.faxis, s.fcode)

"""
    chaintimes(s::SiteLookup, c::SiteFreqChain)
    chainfreq(s::SiteLookup, c::SiteFreqChain)

The time-axis values of a chain's points, and the chain's frequency-channel center.
"""
chaintimes(s::SiteLookup, c::SiteFreqChain) = s.taxis[c.tcodes]
chainfreq(s::SiteLookup, c::SiteFreqChain) = s.faxis[c.fcode]

"""
    sitechains(s::SiteLookup)
    sitechains(s::SiteLookup, site::Symbol)

The per-site chain structure of the lookup: a `NamedTuple` keyed by the *real* site
symbols, whose entries are tuples of [`SiteFreqChain`](@ref Comrade.SiteFreqChain), one
per frequency channel the site observed (ascending in frequency). The two-argument form
returns a single site's tuple.
"""
sitechains(s::SiteLookup) = s.chains
sitechains(s::SiteLookup, site::Symbol) = getproperty(s.chains, site)

"""
    siteindices(s::SiteLookup, site::Symbol)

All flat storage indices belonging to `site`, concatenated across its frequency chains.
"""
siteindices(s::SiteLookup, site::Symbol) =
    reduce(vcat, map(c -> c.inds, sitechains(s, site)))

"""
    npoints(s::SiteLookup)

The total number of parameter points (the length of the flat storage).
"""
npoints(s::SiteLookup) = length(s.tcode)

EnzymeRules.inactive(::typeof(times), ::SiteLookup) = nothing
EnzymeRules.inactive(::typeof(frequencies), ::SiteLookup) = nothing
EnzymeRules.inactive(::typeof(sites), ::SiteLookup) = nothing
EnzymeRules.inactive(::typeof(sitechains), args...) = nothing
EnzymeRules.inactive(::typeof(siteindices), args...) = nothing

function sitemap!(f, out::AbstractArray, gains::AbstractArray, slook::SiteLookup)
    foreach(values(sitechains(slook))) do chs
        foreach(chs) do c
            ysite = @view gains[c.inds]
            outsite = @view out[c.inds]
            outsite .= f.(ysite)
        end
    end
    return out
end

function sitemap(f, gains::AbstractArray{T}, slook::SiteLookup) where {T}
    out = similar(gains)
    sitemap!(f, out, gains, slook)
    return out
end

# Assemble the chain structure from the axes + per-point codes: per (site, channel), the
# point positions sorted ascending in (axis) time.
function _build_sitelookup(taxis, tlo, thi, faxis, flo, fhi, saxis, tcode, fcode, scode)
    chains = map(eachindex(saxis)) do si
        sinds = findall(==(si), scode)
        fcs = sort!(unique(fcode[sinds]))
        Tuple(
            map(fcs) do fc
                finds = sinds[findall(==(fc), fcode[sinds])]
                # ascending in time; `sortperm` is stable so ties keep storage order
                perm = sortperm(tcode[finds])
                inds = finds[perm]
                SiteFreqChain(inds, tcode[inds], fc)
            end
        )
    end
    return SiteLookup(
        NamedTuple{Tuple(saxis)}(Tuple(chains)),
        taxis, tlo, thi, faxis, flo, fhi, saxis, tcode, fcode, scode
    )
end

"""
    SiteArray(arr, sitelookup::SiteLookup)

Construct a site array with the entries `arr` and the site ordering implied by
`sitelookup`.
"""
function SiteArray(a::AbstractArray, map::SiteLookup)
    return SiteArray(a, times(map), frequencies(map), sites(map))
end

function SiteArray(
        data::SiteArray{T, N},
        times::AbstractArray{<:Number, N},
        frequencies::AbstractArray{<:Number, N},
        sites::AbstractArray{<:Any, N}
    ) where {T, N}
    return data
end

function SiteArray(data::SiteArray, ::SiteLookup)
    return data
end
