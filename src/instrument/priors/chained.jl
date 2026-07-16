# The generic chained-prior assembly: everything needed to turn per-(site, frequency)
# chain *distributions* (see `chaindist` in independent.jl) into an observed instrument
# prior over the flat SiteArray storage. This layer owns the index maps between the flat
# storage and the per-chain vectors, the scatter/gather buffers, the merged per-site
# hyperpriors, and the composed transport nodes. It never needs to change when a new
# site prior is added: a site prior only supplies an ordinary `Distributions.Distribution`
# per chain (plus optionally `condition` for exact conditioning and `PT.transport_node`
# for its latent map).

# ----- flat-storage gather/scatter ---------------------------------------------
#
# On the CPU a chain's slice of the flat vector is a (free) view; under Reactant views
# with vector indices do not trace, so ComradeReactantExt overloads these to a
# materialized `stablehlo.gather`/`scatter` (wrapped in `@allowscalar` for the bounds
# check, like `fill_partially_fixed!`).

rgather(x::AbstractVector, inds) = view(x, inds)
function rscatter!(y::AbstractVector, inds, vals)
    y[inds] = vals
    return y
end

# ----- conditioning ordinary distributions ---------------------------------------

"""
    condition(d::Distributions.Distribution, pos, vals)

Exactly condition the multivariate distribution `d` on the values `vals` at the
coordinate positions `pos`, returning a distribution *over the full vector* whose
log-density is the conditional density of the free coordinates (evaluated at a vector
whose fixed positions hold the conditioning values) and whose samples have the fixed
values filled in.

This is how reference-antenna fixing ([`SEFDReference`](@ref), [`SingleReference`](@ref))
and anchoring reach a chained site prior. Implemented for independent products (the fixed
components are dropped) and [`GaussMarkovChain`](@ref Comrade.GaussMarkovChain) (exact
Markov conditioning). A new site prior whose chain distribution can be conditioned in
closed form should add a method; distributions without one simply cannot be combined
with `refant`/anchoring.
"""
function condition(d::Dists.Distribution, pos, vals)
    return throw(
        ArgumentError(
            "$(typeof(d)) does not support exact conditioning, which is required for " *
                "reference-antenna fixing and anchored priors. Define " *
                "`Comrade.condition(::$(nameof(typeof(d))), pos, vals)` or use a site " *
                "prior whose chain distribution supports it (e.g. IIDSitePrior or " *
                "GaussMarkovSitePrior)."
        )
    )
end

function condition(d::Dists.Product, pos::AbstractVector{<:Integer}, vals::AbstractVector{<:Real})
    freepos = setdiff(eachindex(d.v), pos)
    if isempty(freepos)
        v = Vector{Float64}(undef, length(d))
        v[pos] = vals
        return FullyConditionedChain(v)
    end
    # keep the plain product (see `chaindist(::IIDSitePrior, ts)`) — folding overloads
    # of `product_distribution` would lose the per-point transforms
    return PartiallyConditionedDist(
        Dists.Product(d.v[freepos]),
        freepos, collect(Int, pos), collect(float.(vals))
    )
end

# Every point of the chain is conditioned (e.g. an IID override on the reference site):
# the chain carries no prior term and no latent coordinates; the transport just emits
# the conditioning values.
struct FullyConditionedChain <: Dists.ContinuousMultivariateDistribution
    vals::Vector{Float64}
end

Base.length(d::FullyConditionedChain) = length(d.vals)
Base.eltype(::FullyConditionedChain) = Float64
Dists.sampler(d::FullyConditionedChain) = d
Dists.logpdf(d::FullyConditionedChain, x::AbstractVector{<:Number}) = zero(eltype(x))
function Dists._rand!(::Random.AbstractRNG, d::FullyConditionedChain, x::AbstractVector{<:Real})
    x .= d.vals
    return x
end

struct FullyConditionedTransform <: TV.VectorTransform
    vals::Vector{Float64}
end

TV.dimension(::FullyConditionedTransform) = 0
TV.inverse_eltype(::FullyConditionedTransform, x::Type) = eltype(x)

function TV.transform_with(flag::TV.LogJacFlag, t::FullyConditionedTransform, x, index)
    y = similar(x, length(t.vals))
    @inbounds for (i, v) in enumerate(t.vals)
        rsetindex!(y, v, i)
    end
    return y, TV.logjac_zero(flag, eltype(x)), index
end

TV.inverse_at!(x::AbstractArray, index, ::FullyConditionedTransform, y) = index

struct StdFullyConditionedTransform <: PT.AbstractTransport
    vals::Vector{Float64}
end

PT.dimension(::StdFullyConditionedTransform) = 0
PT.pback_eltype(::StdFullyConditionedTransform) = Float64

function PT.pfwd_step(t::StdFullyConditionedTransform, x, index)
    y = similar(x, float(eltype(x)), length(t.vals))
    @inbounds for (i, v) in enumerate(t.vals)
        rsetindex!(y, v, i)
    end
    return y, index
end

PT.pback_step!(x::AbstractVector, index, ::StdFullyConditionedTransform, y) = index

PT.transport_node(d::FullyConditionedChain, ::PT.TVFlat) = FullyConditionedTransform(d.vals)
PT.transport_node(d::FullyConditionedChain, ::PT.AbstractStdDist) = StdFullyConditionedTransform(d.vals)

EnzymeRules.inactive_type(::Type{FullyConditionedChain}) = true
EnzymeRules.inactive_type(::Type{FullyConditionedTransform}) = true
EnzymeRules.inactive_type(::Type{StdFullyConditionedTransform}) = true

# validation hook run at construction so a misconfigured chain (e.g. an unanchored
# Wiener prior) fails at `set_array` time with a clear error instead of at first use
check_chaindist(::Dists.Distribution) = nothing
check_chaindist(d::GaussMarkovChain) = _check_proper(d)

# ----- per-chain log-density / sampling ------------------------------------------
#
# `_cd_logpdf(dist, xc, hp)` is the hyperparameter-aware chain log-density: the generic
# method ignores `hp` (correct for any plain wrapped distribution), GaussMarkovChain
# materializes its process from it. The Product method is a scalar loop so it traces
# under Reactant (the datum loop unrolls on the host; `Distributions`' vectorized product
# logpdf does not trace).

@inline _cd_logpdf(d::Dists.Distribution, x, hp) = Dists.logpdf(d, x)
function _cd_logpdf(d::Dists.Product, x, hp)
    v = d.v
    return sum(k -> Dists.logpdf(v[k], rgetindex(x, k)), eachindex(v); init = zero(eltype(x)))
end
@inline _cd_logpdf(d::PartiallyConditionedDist, x, hp) =
    _cd_logpdf(d.dist, rgather(x, d.variate_index), hp)

# rand is host-only (never traced), so the plain Distributions API suffices
_cd_rand!(rng::Random.AbstractRNG, d::Dists.Distribution, x, hp) = Dists.rand!(rng, d, x)

@inline _selecthp(hp::NamedTuple, ::Val{nothing}) = (;)
@inline _selecthp(hp::NamedTuple, ::Val{K}) where {K} = getproperty(hp, K)

# ----- the assembled distribution over the flat storage --------------------------

# One chain of the assembled prior: the chain distribution, the flat storage indices it
# owns (ascending in time), and the site whose hyperparameters it reads (`Val(nothing)`
# when it has none). Everything here is constant during sampling.
struct ChainSlot{D <: Dists.Distribution, I <: AbstractVector{<:Integer}, K}
    dist::D
    inds::I
    hpsel::Val{K}
end

EnzymeRules.inactive_type(::Type{<:ChainSlot}) = true

"""
    ChainedArrayDist

The observed distribution over the full instrument parameter vector implied by an
`ArrayPrior` containing chained site priors: an independent chain distribution per
(site, frequency channel), each scattered into the flat storage through its index map.
Internal type constructed by `ObservedArrayPrior`.
"""
struct ChainedArrayDist{C <: Tuple} <: Dists.ContinuousMultivariateDistribution
    chains::C
    len::Int
end

Base.length(d::ChainedArrayDist) = d.len
Base.eltype(::ChainedArrayDist) = Float64
Dists.sampler(d::ChainedArrayDist) = d

chainslots(d::ChainedArrayDist) = d.chains
EnzymeRules.inactive(::typeof(chainslots), args...) = nothing

@inline _slot_logpdf(s::ChainSlot, x, hp) =
    _cd_logpdf(s.dist, rgather(x, s.inds), _selecthp(hp, s.hpsel))

function _chain_logpdf(d::ChainedArrayDist, x::AbstractVector, hp::NamedTuple)
    ls = map(s -> _slot_logpdf(s, x, hp), chainslots(d))
    return sum(ls)
end

Dists.logpdf(d::ChainedArrayDist, x::AbstractVector{<:Number}) = _chain_logpdf(d, parent(x), (;))

function _rand_slot!(rng::Random.AbstractRNG, x, s::ChainSlot, hp)
    xc = Vector{Float64}(undef, length(s.dist))
    _cd_rand!(rng, s.dist, xc, _selecthp(hp, s.hpsel))
    rscatter!(x, s.inds, xc)
    return x
end

function _rand_chains!(rng::Random.AbstractRNG, x::AbstractVector, d::ChainedArrayDist, hp::NamedTuple)
    foreach(s -> _rand_slot!(rng, x, s, hp), chainslots(d))
    return x
end

function Dists._rand!(rng::Random.AbstractRNG, d::ChainedArrayDist, x::AbstractVector{<:Real})
    return _rand_chains!(rng, x, d, (;))
end

# ----- the assembled transport nodes ----------------------------------------------
#
# The latent layout is [hyperparameters; chain 1 coords; chain 2 coords; ...]. Each
# chain's transport node is built once from its distribution through the ordinary
# `PT.transport_node` interface — so any wrapped distribution that ProbabilityTransports
# knows gets its transport for free — and evaluated through the hyperparameter-aware
# shims (`_node_transform_with` etc., see markov_chain.jl). Every chain node emits the
# chain's *full* value vector (conditioned points included), so scattering the chains
# covers the whole flat storage.

struct NodeSlot{T, I <: AbstractVector{<:Integer}, K}
    node::T
    inds::I
    hpsel::Val{K}
end

EnzymeRules.inactive_type(::Type{<:NodeSlot}) = true

function _node_slots(d::ChainedArrayDist, space)
    return map(chainslots(d)) do s
        NodeSlot(PT.transport_node(s.dist, space), s.inds, s.hpsel)
    end
end

_slots_flat_dim(slots::Tuple) = reduce(+, map(s -> TV.dimension(s.node), slots); init = 0)
_slots_std_dim(slots::Tuple) = reduce(+, map(s -> PT.dimension(s.node), slots); init = 0)

# thread (logjac, index) over the heterogeneous slot tuple with `Base.tail` recursion
# (TV tuple-transform precedent)
@inline _slots_transform!(flag, y, x, index, ::Tuple{}, hp) = TV.logjac_zero(flag, eltype(x)), index
@inline function _slots_transform!(flag, y, x, index, slots::Tuple, hp)
    s = first(slots)
    yc, ℓ1, index = _node_transform_with(flag, s.node, x, index, _selecthp(hp, s.hpsel))
    rscatter!(y, s.inds, yc)
    ℓ2, index = _slots_transform!(flag, y, x, index, Base.tail(slots), hp)
    return ℓ1 + ℓ2, index
end

@inline _slots_inverse!(x, index, y, ::Tuple{}, hp) = index
@inline function _slots_inverse!(x, index, y, slots::Tuple, hp)
    s = first(slots)
    index = _node_inverse_at!(x, index, s.node, rgather(y, s.inds), _selecthp(hp, s.hpsel))
    return _slots_inverse!(x, index, y, Base.tail(slots), hp)
end

@inline _slots_pfwd!(y, x, index, ::Tuple{}, hp) = index
@inline function _slots_pfwd!(y, x, index, slots::Tuple, hp)
    s = first(slots)
    yc, index = _node_pfwd(s.node, x, index, _selecthp(hp, s.hpsel))
    rscatter!(y, s.inds, yc)
    return _slots_pfwd!(y, x, index, Base.tail(slots), hp)
end

@inline _slots_pback!(x, index, y, ::Tuple{}, hp) = index
@inline function _slots_pback!(x, index, y, slots::Tuple, hp)
    s = first(slots)
    index = _node_pback!(x, index, s.node, rgather(y, s.inds), _selecthp(hp, s.hpsel))
    return _slots_pback!(x, index, y, Base.tail(slots), hp)
end

# --- flat node (fixed hyperparameters); wrapped in `InstrumentTransform` by
# `ObservedArrayPrior`'s generic node, so it returns the plain full-length vector.
struct ChainedTransform{S <: Tuple} <: TV.VectorTransform
    slots::S
    len::Int
    dim::Int
end

TV.dimension(t::ChainedTransform) = t.dim
TV.inverse_eltype(::ChainedTransform, x::Type) = eltype(x)

function TV.transform_with(flag::TV.LogJacFlag, t::ChainedTransform, x, index)
    y = similar(x, t.len)
    ℓ, index = _slots_transform!(flag, y, x, index, t.slots, (;))
    return y, ℓ, index
end

function TV.inverse_at!(x::AbstractArray, index, t::ChainedTransform, y::AbstractVector)
    return _slots_inverse!(x, index, parent(y), t.slots, (;))
end

function PT.transport_node(d::ChainedArrayDist, space::PT.TVFlat)
    slots = _node_slots(d, space)
    return ChainedTransform(slots, d.len, _slots_flat_dim(slots))
end

# --- Std node (fixed hyperparameters); wrapped in `StdInstrumentTransform`.
struct StdChainedTransform{S <: Tuple} <: PT.AbstractTransport
    slots::S
    len::Int
    dim::Int
end

PT.dimension(t::StdChainedTransform) = t.dim
PT.pback_eltype(::StdChainedTransform) = Float64

function PT.pfwd_step(t::StdChainedTransform, x, index)
    y = similar(x, float(eltype(x)), t.len)
    index = _slots_pfwd!(y, x, index, t.slots, (;))
    return y, index
end

function PT.pback_step!(x::AbstractVector, index, t::StdChainedTransform, y)
    return _slots_pback!(x, index, parent(y), t.slots, (;))
end

function PT.transport_node(d::ChainedArrayDist, space::PT.AbstractStdDist)
    slots = _node_slots(d, space)
    return StdChainedTransform(slots, d.len, _slots_std_dim(slots))
end

# ----- the hierarchical observed prior ----------------------------------------

"""
    HierarchicalSample

The shape of a single hierarchical instrument-parameter sample,
`(params = SiteArray, hyperparams = NamedTuple)`, produced by a chained site prior with
fitted hyperparameters (e.g. [`GaussMarkovSitePrior`](@ref)). Downstream code should not
destructure this shape by hand — use [`getparams`](@ref) and [`gethyperparams`](@ref),
which are also defined (as identity/empty) for plain `SiteArray` samples.
"""
const HierarchicalSample = NamedTuple{(:params, :hyperparams)}

"""
    getparams(x)

Extract the instrument parameter values from a single instrument-parameter sample. For
hierarchical samples ([`HierarchicalSample`](@ref Comrade.HierarchicalSample), produced by
chained site priors with fitted hyperparameters) this returns `x.params`; for plain
`SiteArray` samples it is the identity.
"""
@inline getparams(x) = x
@inline getparams(x::HierarchicalSample) = x.params

"""
    gethyperparams(x)

Extract the hyperparameter values from a single instrument-parameter sample: the
`hyperparams` NamedTuple of a [`HierarchicalSample`](@ref Comrade.HierarchicalSample), or
an empty NamedTuple for plain samples.
"""
@inline gethyperparams(x) = (;)
@inline gethyperparams(x::HierarchicalSample) = x.hyperparams

"""
    ObservedHierarchicalArrayPrior

The observed prior produced by an `ArrayPrior` whose chained site priors have fitted
hyperparameters. Samples are NamedTuples `(params = SiteArray, hyperparams = NamedTuple)`
and the log-density is the exact chain density of the parameters given the
hyperparameters plus the hyperprior. Internal type.
"""
struct ObservedHierarchicalArrayPrior{D <: ChainedArrayDist, H <: NamedDist, S <: SiteLookup} <: Dists.ContinuousMultivariateDistribution
    dists::D
    hyperprior::H
    sitemap::S
end

Base.length(d::ObservedHierarchicalArrayPrior) = length(d.dists) + length(d.hyperprior)
Base.eltype(::ObservedHierarchicalArrayPrior) = Float64

function Dists.logpdf(d::ObservedHierarchicalArrayPrior, x::NamedTuple)
    return _chain_logpdf(d.dists, parent(x.params), x.hyperparams) +
        Dists.logpdf(d.hyperprior, x.hyperparams)
end

function Dists.rand(rng::Random.AbstractRNG, d::ObservedHierarchicalArrayPrior)
    hp = rand(rng, d.hyperprior)
    x = Vector{Float64}(undef, length(d.dists))
    _rand_chains!(rng, x, d.dists, hp)
    return (params = SiteArray(x, d.sitemap), hyperparams = hp)
end

function Dists.rand(rng::Random.AbstractRNG, d::ObservedHierarchicalArrayPrior, n::Int)
    return map(_ -> rand(rng, d), 1:n)
end

# Hierarchical transport: the hyperparameter coordinates come *first* in the latent
# vector so the chains can materialize their processes before coloring. The node's shape
# depends only on the sitemap, the conditioned indices, and *which* fields are fitted —
# never on hyperparameter values — preserving the static-transport invariant.
struct WhitenedHierarchicalTransform{H, S <: Tuple, L <: SiteLookup} <: TV.VectorTransform
    hnode::H       # TV node of the hyperprior
    slots::S
    site_map::L
    len::Int
    dim::Int
end

TV.dimension(t::WhitenedHierarchicalTransform) = t.dim

function TV.transform_with(flag::TV.LogJacFlag, t::WhitenedHierarchicalTransform, x, index)
    hp, ℓh, index = TV.transform_with(flag, t.hnode, x, index)
    y = similar(x, t.len)
    ℓc, index = _slots_transform!(flag, y, x, index, t.slots, hp)
    return (params = SiteArray(y, t.site_map), hyperparams = hp), ℓh + ℓc, index
end

function TV.inverse_at!(x::AbstractArray, index, t::WhitenedHierarchicalTransform, y::NamedTuple)
    index = TV.inverse_at!(x, index, t.hnode, y.hyperparams)
    return _slots_inverse!(x, index, parent(y.params), t.slots, y.hyperparams)
end

function TV.inverse_eltype(t::WhitenedHierarchicalTransform, ::Type{T}) where {T <: NamedTuple}
    return promote_type(
        eltype(fieldtype(T, :params)), TV.inverse_eltype(t.hnode, fieldtype(T, :hyperparams))
    )
end

function PT.transport_node(d::ObservedHierarchicalArrayPrior, space::PT.TVFlat)
    hnode = PT.transport_node(d.hyperprior, space)
    slots = _node_slots(d.dists, space)
    return WhitenedHierarchicalTransform(
        hnode, slots, d.sitemap, d.dists.len, TV.dimension(hnode) + _slots_flat_dim(slots)
    )
end

struct StdWhitenedHierarchicalTransform{H, S <: Tuple, L <: SiteLookup, P <: PT.AbstractStdDist} <: PT.AbstractTransport
    hnode::H       # PT node of the hyperprior
    slots::S
    site_map::L
    space::P
    len::Int
    dim::Int
end

PT.dimension(t::StdWhitenedHierarchicalTransform) = t.dim
PT.pback_eltype(::StdWhitenedHierarchicalTransform) = Float64

function PT.pfwd_step(t::StdWhitenedHierarchicalTransform, x, index)
    hp, index = PT.pfwd_step(t.hnode, x, index)
    y = similar(x, float(eltype(x)), t.len)
    index = _slots_pfwd!(y, x, index, t.slots, hp)
    return (params = SiteArray(y, t.site_map), hyperparams = hp), index
end

function PT.pback_step!(x::AbstractVector, index, t::StdWhitenedHierarchicalTransform, y::NamedTuple)
    index = PT.pback_step!(x, index, t.hnode, y.hyperparams)
    return _slots_pback!(x, index, parent(y.params), t.slots, y.hyperparams)
end

function PT.transport_node(d::ObservedHierarchicalArrayPrior, space::PT.AbstractStdDist)
    hnode = PT.transport_node(d.hyperprior, space)
    slots = _node_slots(d.dists, space)
    return StdWhitenedHierarchicalTransform(
        hnode, slots, d.sitemap, space, d.dists.len, PT.dimension(hnode) + _slots_std_dim(slots)
    )
end

# ----- construction from ArrayPrior --------------------------------------------

function build_chain_observed(d::ArrayPrior, site_dists::NamedTuple, smap::SiteLookup, array)
    d.centroid_station === nothing || throw(
        ArgumentError(
            "centroid_station is only supported when every site prior is a product " *
                "(e.g. IIDSitePrior); it cannot be combined with chained site priors."
        )
    )

    finds, vals = reference_indices(array, smap, d.refant)
    refvals = vals === nothing ? Float64[] : collect(float.(vals))
    refmap = Dict{Int, Float64}(zip(collect(Int, finds), refvals))

    # Each (site, frequency) chain contributes one slot: its chain distribution from the
    # site-prior interface (`chaindist`), exactly conditioned on any reference-fixed
    # points plus the prior's own anchor points (reference fixing wins duplicates).
    slots = Any[]
    for (s, chs) in pairs(sitechains(smap))
        sp = getproperty(site_dists, s)
        hpk = isempty(hyperprior(sp)) ? nothing : s
        for c in chs
            inds = collect(Int, c.inds)
            pos = Int[]
            pvals = Float64[]
            for (j, i) in enumerate(inds)
                if haskey(refmap, i)
                    push!(pos, j)
                    push!(pvals, refmap[i])
                end
            end
            ainds, avals = anchor_points(sp, c)
            for (i, v) in zip(ainds, avals)
                j = findfirst(==(i), inds)
                (j === nothing || j ∈ pos) && continue
                push!(pos, j)
                push!(pvals, v)
            end
            dist = chaindist(sp, chaintimes(smap, c))
            length(dist) == length(inds) || throw(
                ArgumentError(
                    "chaindist for site $s returned a distribution of length " *
                        "$(length(dist)); expected $(length(inds)) (one per chain point)."
                )
            )
            isempty(pos) || (dist = condition(dist, pos, pvals))
            check_chaindist(dist)
            push!(slots, ChainSlot(dist, inds, Val(hpk)))
        end
    end
    chains = Tuple(slots)

    # Assemble the merged hyperprior: every site gets its own hyperparameters, nested
    # under the site name. The site prior in the ArrayPrior default is a per-site
    # template, so its fitted fields are replicated for each site that uses it.
    # Multi-frequency chains of the same site share that site's hyperparameters.
    hp = (;)
    for (s, sp) in pairs(site_dists)
        h = hyperprior(sp)
        isempty(h) && continue
        hp = merge(hp, NamedTuple{(s,)}((h,)))
    end

    dists = ChainedArrayDist(chains, npoints(smap))
    isempty(hp) && return ObservedArrayPrior(dists, smap)
    return ObservedHierarchicalArrayPrior(dists, NamedDist(hp), smap)
end
