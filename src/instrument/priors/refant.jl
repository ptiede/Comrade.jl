abstract type ReferencingScheme end

export NoReference, SingleReference, SEFDReference

struct NoReference <: ReferencingScheme end

"""
    SingleReference(site::Symbol, val)

Selects a single reference site for all scans. The value of the site is set to `val`.
"""
struct SingleReference{T} <: ReferencingScheme
    site::Symbol
    value::T
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

reference_indices(::AbstractArrayConfiguration, st::SiteLookup, ::NoReference) = [], nothing
function reference_indices(::AbstractArrayConfiguration, st::SiteLookup, p::SingleReference)
    inds = findall(==(p.site), st.sites)
    return inds, Fill(p.value, length(inds))
end

function reference_indices(array::AbstractArrayConfiguration, st::SiteLookup, r::SEFDReference)
    tarr = array.tarr
    t = unique(st.times)
    sefd = NamedTuple{Tuple(tarr.sites)}(Tuple(tarr.SEFD1 .+ tarr.SEFD2))
    fixedinds = map(eachindex(t)) do i
        inds = findall(==(t[i]), st.times)
        sites = Tuple(st.sites[inds])
        @assert length(sites) < length(sefd) "Error in reference site generation. Too many sites"
        sp = select(sefd, sites)
        _, ind = findmin(values(sp))
        return inds[ind]
    end
    return fixedinds, Fill(r.value, length(fixedinds))
end
