abstract type ReferencingScheme end

export NoReference, SingleReference, RandomReference, SEFDReference

struct NoReference <: ReferencingScheme end

struct SingleReference{T} <: ReferencingScheme
    site::Symbol
    value::T
end

"""
    SingleReference(site::Symbol, val::Number)

Use a single site as a reference. The sites gain will be set equal to `val`.
"""
SingleReference(sites::Symbol, scheme::Number) = SingleReference(sites, FixedSeg(scheme))

struct RandomReference{T} <: ReferencingScheme
    value::T
end

"""
    RandomReference(val::Number)

For each timestamp select a random reference sites whose sites gain will be set to `val`.

## Notes
This is useful when there isn't a single site available for all scans and you want to split
up the choice of reference site. We recommend only using this option for Stokes I fitting.
"""
RandomReference(val::Number) = RandomReference(FixedSeg(val))

struct SEFDReference{T} <: ReferencingScheme
    value::T
    offset::Int
end

"""
    SiteOrderReference(val::Number, sefd_index = 1)

Selects the reference site based on the SEFD of each telescope, where the smallest SEFD
is preferentially selected. The reference gain is set to `val` and the user can select to
use the `n` lowest SEFD site by passing `sefd_index = n`.

## Notes
This is done on a per-scan basis so if a site is missing from a scan the next highest SEFD
site will be used.
"""
SEFDReference(val::Number, offset::Int=0) = SEFDReference{typeof(FixedSeg(val))}(FixedSeg(val), offset)

reference_sites(st::ScanTable, ::NoReference) = Fill(NoReference(), length(st))
reference_sites(st::ScanTable, p::SingleReference) = Fill(p, length(st))

function reference_sites(st::ScanTable, p::RandomReference)
    return map(1:length(st)) do i
        s = st[i]
        ref = rand(sites(s))
        return SingleReference(ref, p.scheme)
    end
end

function reference_sites(st::ScanTable, p::SEFDReference)
    tarr = st.obs.config.tarr
    sefd = NamedTuple{Tuple(tarr.sites)}(Tuple(tarr.SEFD1 .+ tarr.SEFD2))
    map(1:length(st)) do i
        s = st[i]
        sites = sites(s)
        sp = select(sefd, sites)
        ind = findmin(values(sp))[2]
        indo = (ind+p.offset > length(sites)) ? (ind+p.offset)%length(sites) : ind+p.offset
        return SingleReference(sites[indo], p.scheme)
    end
end
