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
SEFDReference(val::Number, offset::Int=0) = SEFDReference{typeof(FixedSeg(val))}(FixedSeg(val), offset)

reference_sites(::AbstractArrayConfiguration, st::SiteLookup, ::NoReference) = Fill(NoReference(), length(st))
reference_sites(::AbstractArrayConfiguration, st::SiteLookup, p::SingleReference) = Fill(p, length(st))

function reference_sites(array::AbstractArrayConfiguration, st::SiteLookup, p::SEFDReference)
    tarr = array.tarr
    t = unique(sitemap.times)
    sefd = NamedTuple{Tuple(tarr.sites)}(Tuple(tarr.SEFD1 .+ tarr.SEFD2))
    map(eachindex(t)) do i
        inds = findall(==(t[i]), st.times)
        sites = Tuple(sitmap.sites[inds])
        @assert length(sites) < length(sefd) "Error in reference site generation. Too many sites"
        sp = select(sefd, sites)
        _, ind = findmin(values(sp))[2]
        return SingleReference(sites[ind], p.value), inds[ind]
    end
end

struct PartiallyConditionedDist{D<:Distributions.ContinuousMultivariateDistribution, I, F} <: Distributions.ContinuousMultivariateDistribution
    dist::D
    variate_index::I
    fixed_index::I
    fixed_values::F
end

Base.length(d::PartiallyConditionedDist) = length(d.variate_index) + length(d.fix)
Base.size(d::PartiallyConditionedDist) = size(d.variate_index) .+ size(d.fixed_index)
Base.eltype(d::PartiallyConditionedDist) = eltype(d.dist)

function Distributions._logpdf(d::PartiallyConditionedDist, x)
    xv = @view x[d.variate_index]
    return logpdf(d.dist, xv)
end

function Distributions._rand!(rng::AbstractRNG, d::PartiallyConditionedDist, x::AbstractArray)
    rand!(rng, d.dist, @view(x[d.variate_index]))
    # Now adjust the other indices
    x[d.fixed_index] .= d.fixed_values
    return x
end
