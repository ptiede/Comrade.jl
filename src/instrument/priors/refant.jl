abstract type ReferencingScheme end

export NoReference, SingleReference, MultiReference, SEFDReference

struct NoReference <: ReferencingScheme end

"""
    SingleReference(site::Symbol, val)

Selects a single reference site for all scans. The value of the site is set to `val`.
"""
struct SingleReference{T} <: ReferencingScheme
    site::Symbol
    value::T
end


"""
    MultiReference(sites, val)

Selects several reference sites at once, fixing every entry of each site in `sites` to
`val`. `sites` is any iterable of site names.

One pin removes one gauge freedom, so this is the scheme for a parameter that carries
more than one. A track-segmented phase offset is the usual case: the array partitions
into groups whose absolute phase level the schedule never links, and each group needs its
own pin. That happens whenever a subset of the array observes for a stretch with no
station in common with the site that is already pinned.

## Notes
This pins exactly the sites it is given and does not check that the result leaves the
model identified. Too few pins leave a flat direction in the posterior; too many
over-constrain a scan, since a scan carries only one global phase.
"""
struct MultiReference{S <: AbstractVector{Symbol}, T} <: ReferencingScheme
    sites::S
    value::T
    function MultiReference{S, T}(sites, value) where {S, T}
        ss = convert(S, sites)::S
        isempty(ss) && throw(ArgumentError("MultiReference needs at least one site"))
        allunique(ss) ||
            throw(ArgumentError("MultiReference sites must be distinct, got $ss"))
        return new{S, T}(ss, value)
    end
end

MultiReference{S}(sites, value) where {S} = MultiReference{S, typeof(value)}(sites, value)
MultiReference(sites::AbstractVector{Symbol}, value) =
    MultiReference{typeof(sites), typeof(value)}(sites, value)
MultiReference(sites, value) = MultiReference(Symbol[Symbol(s) for s in sites], value)

#    EntryReference(indices, value)
#
# Pins the parameter entries at `indices` — positions in the parameter's `SiteLookup` — to
# `value`. Unlike the site-named schemes this one carries no notion of which site or stamp
# it pins, so it is meaningful only alongside the `SiteLookup` the indices were computed
# from. `gauge_pins` returns such indices and `CompositeReference` is how they are added to
# the scheme a user already asked for. Internal type.
struct EntryReference{I <: AbstractVector{Int}, T} <: ReferencingScheme
    indices::I
    value::T
    function EntryReference{I, T}(indices, value) where {I, T}
        ii = convert(I, indices)::I
        allunique(ii) ||
            throw(ArgumentError("EntryReference indices must be distinct, got $ii"))
        return new{I, T}(ii, value)
    end
end

EntryReference{I}(indices, value) where {I} = EntryReference{I, typeof(value)}(indices, value)
EntryReference(indices::AbstractVector{Int}, value) =
    EntryReference{typeof(indices), typeof(value)}(indices, value)
EntryReference(indices, value) = EntryReference(collect(Int, indices), value)

#    CompositeReference(schemes...)
#
# Applies every scheme in `schemes` to the same parameter and pins the union of their
# indices. An index two schemes both claim must be given the same value by both. Internal
# type.
struct CompositeReference{S <: Tuple} <: ReferencingScheme
    schemes::S
end

CompositeReference(schemes::ReferencingScheme...) = CompositeReference(schemes)

# The value a scheme pins its entries to; `nothing` when it pins nothing.
reference_value(::NoReference) = nothing
reference_value(r::ReferencingScheme) = r.value
function reference_value(r::CompositeReference)
    for s in r.schemes
        v = reference_value(s)
        v === nothing || return v
    end
    return nothing
end

struct SEFDReference{T} <: ReferencingScheme
    value::T
    offset::Int
end

"""
    SEFDReference(val::Number, sefd_index = 1)

Selects the reference site based on the SEFD of each telescope, where the smallest SEFD
is preferentially selected. The reference gain is set to `val` and the user can select to
use the `n` lowest SEFD site by passing `sefd_index = n`.

## Notes
This is done on a per-scan basis so if a site is missing from a scan the next highest SEFD
site will be used.
"""
SEFDReference(val::Number) = SEFDReference(val, 0)

# Only SEFDReference reads the array configuration; the schemes that name their sites
# outright ignore it.
reference_indices(_, st::SiteLookup, ::NoReference) = [], nothing
function reference_indices(_, st::SiteLookup, p::SingleReference)
    inds = findall(==(p.site), st.sites)
    return inds, fill(p.value, length(inds))
end

function reference_indices(_, st::SiteLookup, p::MultiReference)
    inds = Int[]
    for s in p.sites
        si = findall(==(s), st.sites)
        # A name absent from this parameter's site list pins nothing, which would leave
        # the gauge the scheme was asked to fix silently unfixed.
        isempty(si) && throw(
            ArgumentError(
                "MultiReference site $s is not among this parameter's sites " *
                    "$(sort(unique(st.sites)))"
            )
        )
        append!(inds, si)
    end
    sort!(inds)
    return inds, fill(p.value, length(inds))
end

function reference_indices(_, st::SiteLookup, r::EntryReference)
    inds = sort(r.indices)
    for i in inds
        checkindex(Bool, eachindex(st.sites), i) || throw(
            ArgumentError(
                "EntryReference index $i is outside this parameter's " *
                    "$(length(st.sites)) entries"
            )
        )
    end
    return inds, fill(r.value, length(inds))
end

function reference_indices(array, st::SiteLookup, r::CompositeReference)
    parts = [reference_indices(array, st, s) for s in r.schemes]
    filter!(p -> !isempty(first(p)), parts)
    isempty(parts) && return Int[], nothing
    T = mapreduce(p -> eltype(last(p)), promote_type, parts)
    inds = Int[]
    vals = T[]
    for (is, vs) in parts, (i, v) in zip(is, vs)
        j = findfirst(==(i), inds)
        if j === nothing
            push!(inds, i)
            push!(vals, v)
        elseif vals[j] != v
            throw(
                ArgumentError(
                    "CompositeReference pins entry $i (site $(st.sites[i])) to both " *
                        "$(vals[j]) and $v"
                )
            )
        end
    end
    p = sortperm(inds)
    return inds[p], vals[p]
end

function reference_indices(array::AbstractArrayConfiguration, st::SiteLookup, r::SEFDReference)
    tarr = array.tarr
    t = unique(st.times)
    f = unique(st.frequencies)
    sefd = NamedTuple{Tuple(tarr.sites)}(Tuple(tarr.SEFD1 .+ tarr.SEFD2))
    fixedinds = Int[]
    for i in eachindex(t), j in eachindex(f)
        inds = findall(x -> ((st.times[x] == t[i])&&(st.frequencies[x] == f[j])), eachindex(st.times))
        if isempty(inds)
            continue
        end
        sites = Tuple(st.sites[inds])
        @assert length(sites) <= length(sefd) "Error in reference site generation. Too many sites"
        _, ind = findmin(map(s -> sefd[s], sites))
        push!(fixedinds, inds[ind])
    end
    return fixedinds, fill(r.value, length(fixedinds))
end
