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

Selects several reference sites, fixing every entry of each site in `sites` to `val`.
`sites` is any iterable of site names.

Use it for a parameter with more than one gauge freedom. The usual case is a
track-segmented phase offset when the schedule splits the array into groups that share no
station, for example when a subset of the array observes for a stretch without the
reference site. Each such group needs its own reference.

## Notes
This pins the sites it is given and does not check that the result leaves the model
identified. Too few pins leave a flat direction in the posterior; too many
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
# Pins the parameter entries at `indices`, positions in the parameter's `SiteLookup`, to
# `value`. The indices are meaningful only for the `SiteLookup` they were computed from.
# Internal type.
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

# Only SEFDReference reads the array configuration.
reference_indices(_, st::SiteLookup, ::NoReference) = [], nothing
function reference_indices(_, st::SiteLookup, p::SingleReference)
    inds = findall(==(p.site), st.sites)
    return inds, fill(p.value, length(inds))
end

function reference_indices(_, st::SiteLookup, p::MultiReference)
    inds = Int[]
    for s in p.sites
        si = findall(==(s), st.sites)
        # An absent site would pin nothing and leave the gauge unfixed.
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
