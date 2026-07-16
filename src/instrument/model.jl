export InstrumentModel

"""
    $(TYPEDEF)

The abstract instrument model. For a concrete implementation see [`IdealInstrumentModel`](@ref) and
[`InstrumentModel`](@ref).

Any subtype must implement the following methods

 - `set_array(m::AbstractInstrumentModel, array::AbstractArrayConfiguration)`: Sets the array configuration
    for the instrument model `m` and returns the observed instrument model and prior.
 - `apply_instrument(vis, m::AbstractInstrumentModel, x)`: Applies the instrument model `m` to the visibilities
    `vis` given the model parameters `x`.
"""
abstract type AbstractInstrumentModel end

"""
    IdealInstrument(array::AbstractArrayConfiguration)

Constructs an ideal instrument that has no corruptions or feed rotation.
"""
struct IdealInstrumentModel <: AbstractInstrumentModel end

Base.show(io::IO, mime::MIME"text/plain", m::IdealInstrumentModel) = printstyled(io, "IdealInstrumentModel"; color = :light_cyan, bold = true)

@inline apply_instrument(vis, ::IdealInstrumentModel, x) = vis


struct InstrumentModel{J <: AbstractJonesMatrix, PI, P <: PolBasis} <: AbstractInstrumentModel
    jones::J
    prior::PI
    refbasis::P
end

function Base.show(io::IO, ::MIME"text/plain", m::InstrumentModel)
    printstyled(io, "InstrumentModel"; bold = true, color = :light_cyan)
    println(io)
    T = typeof(m.jones)
    ST = split(split(" $T", '{')[1], ".")[end]
    println(io, "  with Jones: ", ST)
    return print(io, "  with reference basis: ", m.refbasis)
end


struct ObservedInstrumentModel{I <: AbstractJonesMatrix, PB <: PolBasis, B, M <: NamedTuple} <: AbstractInstrumentModel
    """
    The abstract instrument model
    """
    instrument::I
    """
    reference basis used to define the ideal basis

    """
    refbasis::PB
    """
    The baseline site lookup for the instrument model
    """
    bsitelookup::B
    """
    Per-visibility metadata (Ti, Fr, and per-slot site vectors) used to build each `JonesPoint`
    """
    vismeta::M
end

# Site lookup is const so we add a method so we can signal
# to Enzyme that it is not differentiable.
sitelookup(x::ObservedInstrumentModel) = x.bsitelookup
instrument(x::ObservedInstrumentModel) = x.instrument
refbasis(x::ObservedInstrumentModel) = x.refbasis
vismeta(x::ObservedInstrumentModel) = x.vismeta
EnzymeRules.inactive(::typeof(sitelookup), args...) = nothing
EnzymeRules.inactive(::typeof(instrument), args...) = nothing
EnzymeRules.inactive(::typeof(refbasis), args...) = nothing
EnzymeRules.inactive(::typeof(vismeta), args...) = nothing

function Base.show(io::IO, mime::MIME"text/plain", m::ObservedInstrumentModel)
    printstyled(io, "ObservedInstrumentModel"; bold = true, color = :light_cyan)
    println(io)
    T = typeof(m.instrument)
    ST = split(split(" $T", '{')[1], ".")[end]
    println(io, "  with Jones: ", ST)
    return print(io, "  with reference basis: ", m.refbasis)
end


"""
    InstrumentModel(jones, prior; refbasis = CirBasis())

Builds an instrument model using the jones matrix `jones`, with priors `prior`.
The reference basis is `refbasis` and is used to define what
the ideal basis is. Namely, the basis that you have the ideal visibilties to be represented in.
For classical VLBI `refbasis = CirBasis` is a good default option, sinc the majority of the
array uses circular feeds. For linear feed arrays like VGOS a user should switch to `LinBasis`,
although failure to do so will not cause any errors, and is just a less efficient representation of the
visibilities.

# Arguments

 - `jones` : The jones matrix that represents the instrument. This is a function that takes in the
    parameters of the instrument and returns a jones matrix. See [`SingleStokesGain`](@ref)
    for a Stokes I example and [`JonesG`](@ref) or [`JonesD`](@ref) for polarized examples.
 - `prior`: A named tuple of [`ArrayPrior`](@ref) that specify what the priors are for each
    component used to construct the jones matrix using the function `jones`


# Optional Arguments
  - `refbasis`: The reference basis used for the computation. The default is `CirBasis()` which are circular feeds.


# Example

A Stokes I example is
```julia-repl
julia> G = SingleStokesGain((x, p)->exp(complex(x.lg, x.pg)))
julia> intprior = (lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            pg = ArrayPrior(IIDSitePrior(ScanSeg(), DiagVonMises(0.0, inv(π^2))))
            )

julia> intm = InstrumentModel(G, intprior)
```

A standard polarized example is
```julia-repl
julia> G = JonesG() do x, p
        gR = exp.(complex(x.lgr, x.gpr))
        gL = gr*exp.(complex(x.lgrat, x.gprat))
        return gR, gL
    end
julia> D = JonesD() do x, p
        dR = complex.(x.dRre, x.dRim)
        dL = complex.(x.dLre, x.dLim)
        return gR, gL
    end
julia> R = JonesR()
julia> J = JonesSandwich(G, D, R)
julia> intprior = (lgr = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)),
                    gpr = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))),
                    lgrat = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)),
                    gprat = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))),
                    dRre = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1)),
                    dRim = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1)),
                    dLre = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1)),
                    dLim = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1))
                    )
julia> intm = InstrumentModel(J, intprior)
```

which construct the gain matrix from R and ratios, and D is the small leakage matrix. [`JonesR`](@ref)
is the *response matrix* that controls how the site responds to the ideal visibility in the reference
basis.
"""
# Every param_map takes (x, p::JonesPoint) as of Comrade 0.12. A one-argument function
# would only fail deep inside the likelihood with a MethodError, so check at construction
# time and point at the migration. Composite/parameter-free types dispatch; the fallback
# checks any `param_map` field of user-defined Jones types.
_check_param_map_arity(j::JonesSandwich) = foreach(_check_param_map_arity, j.matrices)
_check_param_map_arity(::Union{JonesF, JonesR, JonesConst}) = nothing
function _check_param_map_arity(j::AbstractJonesMatrix)
    hasfield(typeof(j), :param_map) || return nothing
    f = getfield(j, :param_map)
    f isa Function || return nothing
    ms = methods(f)
    isempty(ms) && return nothing
    # nargs counts the function itself, so a two-argument method has nargs == 3
    any(m -> m.nargs == 3 || m.isva, ms) && return nothing
    throw(
        ArgumentError(
            "the `param_map` of $(nameof(typeof(j))) must take two arguments `(x, p)` as of " *
                "Comrade 0.12: `x` holds the sampled parameter values and `p::JonesPoint` the " *
                "gain point's metadata (`p.Ti`, `p.Fr`, `p.site`). Update e.g. " *
                "`SingleStokesGain(x -> exp(x.lg))` to `SingleStokesGain((x, p) -> exp(x.lg))`."
        )
    )
end

function InstrumentModel(jones::AbstractJonesMatrix, prior::NamedTuple{N}; refbasis = CirBasis()) where {N}
    _check_param_map_arity(jones)
    return InstrumentModel(jones, prior, refbasis)
end

function InstrumentModel(jones::JonesR; refbasis = CirBasis())
    return InstrumentModel(jones, NamedTuple(), refbasis)
end

function set_array(int::InstrumentModel, array::AbstractArrayConfiguration)
    (; jones, prior, refbasis) = int
    # 1. preallocate and jones matrices
    Jpre = preallocate_jones(jones, array, refbasis)
    # 2. construct the prior with the array you have
    obs = map(x -> ObservedArrayPrior(x, array), prior)
    prior_obs = NamedDist(obs)
    # 3. construct the baseline site map for each prior from its sitemap
    bsitemaps = map(d -> _construct_baselinemap(array, d.sitemap), obs)
    # 4. per-visibility metadata for the JonesPoints (plain vectors: zero lookup cost;
    #    the site columns are indexed by antenna slot)
    bl = array[:sites]
    meta = (
        Ti = float.(array[:Ti]),
        Fr = float.(array[:Fr]),
        s = (first.(bl), last.(bl)),
    )
    intobs = ObservedInstrumentModel(Jpre, refbasis, bsitemaps, meta)
    return intobs, prior_obs
end

function set_array(int::IdealInstrumentModel, ::AbstractArrayConfiguration)
    return (int, NamedDist((;)))
end

struct BaselineSiteLookup{V <: AbstractArray}
    indices_1::V
    indices_2::V
end

function _construct_baselinemap(array::EHTArrayConfiguration, smap::SiteLookup)
    T = array[:Ti]
    F = array[:Fr]
    bl = array[:sites]
    return _construct_baselinemap(T, F, bl, smap)
end

# The position of the (unique) point of chain `c` whose half-open time segment
# [tlo, thi) contains `t`, plus the number of matching points (overlapping segments are
# a data error). The chain's segments are ascending, so a binary search over the segment
# ends replaces the old per-datum linear scan.
@inline function _chain_time_index(smap::SiteLookup, c::SiteFreqChain, t)
    codes = c.tcodes
    ends = CodedAxis(smap.thi, codes)
    j = searchsortedfirst(ends, t)  # first point whose segment end is ≥ t
    idx = 0
    nmatch = 0
    k = j
    @inbounds while k ≤ lastindex(codes) && smap.tlo[codes[k]] ≤ t
        if t < smap.thi[codes[k]]
            nmatch += 1
            idx == 0 && (idx = k)
        end
        k += 1
    end
    return idx, nmatch
end

function _sitefreqindex(smap::SiteLookup, s, t, f)
    haskey(sitechains(smap), s) ||
        throw(AssertionError("$t, $f, $((s)) not found in SiteArray"))
    idx = 0
    nmatch = 0
    for c in sitechains(smap, s)
        (smap.flo[c.fcode] ≤ f < smap.fhi[c.fcode]) || continue
        j, n = _chain_time_index(smap, c, t)
        n == 0 && continue
        idx = c.inds[j]
        nmatch += n
    end
    nmatch > 1 && throw(AssertionError("Multiple indices found for $t, $((s)) in SiteArray"))
    nmatch == 0 && throw(AssertionError("$t, $f, $((s)) not found in SiteArray"))
    return idx
end

function _construct_baselinemap(T, F, bl, smap::SiteLookup)
    ind1 = similar(T, Int)
    ind2 = similar(T, Int)
    for i in eachindex(T, F, bl, ind1, ind2)
        t = T[i]
        f = F[i]
        s1, s2 = bl[i]
        ind1[i] = _sitefreqindex(smap, s1, t, f)
        ind2[i] = _sitefreqindex(smap, s2, t, f)
    end
    return BaselineSiteLookup(ind1, ind2)
end


@inline intout(vis::AbstractArray{<:StokesParams{T}}) where {T} = similar(vis, SMatrix{2, 2, complex(T), 4})
@inline intout(vis::AbstractArray{T}) where {T} = similar(vis, complex(T))
@inline intout(vis::AbstractArray{<:CoherencyMatrix{A, B, T}}) where {A, B, T} = similar(vis, SMatrix{2, 2, complex(T), 4})
@inline intout(vis::StructArray{<:StokesParams{T}}) where {T} = StructArray{SMatrix{2, 2, complex(T), 4}}((similar(vis.I), similar(vis.Q), similar(vis.U), similar(vis.V)))

@inline function apply_instrument(vis, J::ObservedInstrumentModel, x)
    vout = parent(intout(vis))
    # Grab parent arrary so we don't trace through SiteArray since that stuff is constant.
    # `getparams` strips the hyperparameters from hierarchical samples first.
    xint = map(parent ∘ getparams, x.instrument)
    _apply_instrument!(vout, parent(vis), J, xint)
    return vout
end


@inline function _apply_instrument!(vout::AbstractArray, vis, J::ObservedInstrumentModel, xint)
    @inbounds for i in eachindex(vout, vis)
        vout[i] = @inline apply_jones(vis[i], i, J, xint)
    end
    return
end

EnzymeRules.inactive_type(::Type{<:ObservedInstrumentModel}) = true


@inline function apply_instrument(vis, J::ObservedInstrumentModel{<:JonesConst}, x)
    vout = intout(parent(vis))
    vout .= apply_jones.(vis, eachindex(vis), Ref(J), Ref((;)))
    return vout
end

#EnzymeRules.inactive(::typeof(Base.Ref), ::ObservedInstrumentModel) = nothing

# @inline function _apply_instrument!(vout, vis, J::ObservedInstrumentModel, xint)
#     # @inbounds for i in eachindex(vout, vis)
#     #     v = apply_jones(vis[i], i, J, xint)
#     #     vout[i] = v
#     # end
#     vout .= apply_jones.(vis, eachindex(vis), Ref(J), Ref(xint))
#     return nothing
# end

@inline get_indices(bsitemaps, index, ::Val{1}) = map(Base.Fix2(rgetindex, index) ∘ Base.Fix2(getproperty, :indices_1), bsitemaps)
@inline get_indices(bsitemaps, index, ::Val{2}) = map(Base.Fix2(rgetindex, index) ∘ Base.Fix2(getproperty, :indices_2), bsitemaps)
@inline get_params(x::NamedTuple{N}, indices::NamedTuple{N}) where {N} = NamedTuple{N}(map(rgetindex, values(x), values(indices)))
# @inline get_params(x::NamedTuple{N}, indices::NamedTuple{N}) where {N} = NamedTuple{N}(ntuple(i->getindex(x[i], indices[i]), Val(length(N))))

# We need this because Enzyme seems to crash when generating code for this
# TODO try to find MWE and post to Enzyme.jl
EnzymeRules.inactive(::typeof(get_indices), args...) = nothing

@inline function unrollmap(x::NTuple{N}, inds::NTuple{N}) where {N}
    return ntuple(Val(N)) do k
        Base.@_inline_meta
        @inbounds getindex(x[k], inds[k])
    end
end

# Build the JonesPoint of one antenna slot of datum `i`. Overloaded in the Reactant
# extension for traced indices, where the site Symbol cannot be gathered and is `nothing`.
Base.@propagate_inbounds _sitepoint(::Val{N}, meta, i::Integer) where {N} =
    JonesPoint{N}(meta.Ti[i], meta.Fr[i], meta.s[N][i], i)

# The Jones matrix of antenna slot `S` (1 = lhs, 2 = rhs station) of datum `index`:
# gather each parameter's value through the exact per-datum index maps built at
# `set_array` time. Shared by the likelihood hot path (`apply_jones`) and the
# diagnostic forward model (`forward_jones(::ObservedInstrumentModel, xs)`), so what
# you inspect is what corrupts the visibilities.
Base.@propagate_inbounds function slot_jones(J::ObservedInstrumentModel, x::NamedTuple{N}, ::Val{S}, index) where {N, S}
    indices = get_indices(sitelookup(J), index, Val(S))
    params = NamedTuple{N}(unrollmap(values(x), values(indices)))
    p = _sitepoint(Val(S), vismeta(J), index)
    return jonesmatrix(instrument(J), params, p)
end

Base.@propagate_inbounds function apply_jones(v, index, J::ObservedInstrumentModel, x::NamedTuple)
    j1 = slot_jones(J, x, Val(1), index)
    j2 = slot_jones(J, x, Val(2), index)
    vout = _apply_jones(v, j1, j2, refbasis(J))
    return vout
end

"""
    forward_jones(J::ObservedInstrumentModel, xs::NamedTuple)

Construct the Jones matrices of the observed gain model `J` at every point of the densest
element of `xs`, associating the parameters through the exact per-datum index maps built
at `set_array` time — the same lookups the likelihood uses, so mixed segmentations need
no matching. Each output point is evaluated at a representative datum it contains, which
is well defined whenever the coarser parameters are constant across the point's segment
(the usual nested segmentations). This is what [`instrumentmodel(post, θ)`](@ref) calls.
"""
function forward_jones(J::ObservedInstrumentModel, xs::NamedTuple)
    ib = argmax(map(length ∘ times, values(xs)))
    xb = values(xs)[ib]
    bmap = values(sitelookup(J))[ib]
    # representative (datum, antenna slot) of every output point; every point of the
    # lookup's grid contains at least one datum by construction
    repi = zeros(Int, length(xb))
    reps = zeros(Int, length(xb))
    for i in eachindex(bmap.indices_1, bmap.indices_2)
        j1 = bmap.indices_1[i]
        if repi[j1] == 0
            repi[j1] = i
            reps[j1] = 1
        end
        j2 = bmap.indices_2[i]
        if repi[j2] == 0
            repi[j2] = i
            reps[j2] = 2
        end
    end
    vs = [_rep_jones(J, xs, repi[j], reps[j]) for j in eachindex(repi)]
    return SiteArray(vs, times(xb), frequencies(xb), sites(xb))
end

function _rep_jones(J::ObservedInstrumentModel, xs::NamedTuple, i, s)
    s == 1 && return slot_jones(J, xs, Val(1), i)
    return slot_jones(J, xs, Val(2), i)
end


@inline _apply_jones(v::Number, j1, j2, ::B) where {B} = j1 * v * conj(j2)
@inline _apply_jones(v::CoherencyMatrix, j1, j2, ::B) where {B} = j1 * CoherencyMatrix{B, B}(v) * j2'
@inline _apply_jones(v::StokesParams, j1, j2, ::B) where {B} = j1 * CoherencyMatrix{B, B}(v) * j2'


# function ChainRulesCore.rrule(::typeof(apply_instrument), vis, J::ObservedInstrumentModel, x)
#     out = apply_instrument(vis, J, x)
#     px = ProjectTo(x)
#     function _apply_instrument_pb(Δ)
#         bvis = baseimage(vis)
#         bout = baseimage(out)
#         Δout = similar(bout)
#         Δout .= unthunk(Δ)
#         xi = x.instrument
#         dx = ntzero(xi)
#         dvis = zero(bvis)
#         autodiff(Reverse, _apply_instrument!, Const, Duplicated(bout, Δout), Duplicated(bvis, dvis), Const(J), Duplicated(xi, dx))
#         return NoTangent(), UnstructuredMap(dvis, axisdims(vis)), NoTangent(), px((;instrument = dx))
#     end
#     return out, _apply_instrument_pb
# end

# function ChainRulesCore.rrule(::typeof(apply_instrument), vis, J::ObservedInstrumentModel{<:Union{JonesR, JonesF}}, x)
#     out = apply_instrument(vis, J, x)
#     function _apply_instrument_pb(Δ)
#         Δout = similar(out)
#         Δout .= unthunk(Δ)
#         dvis = zero(vis)
#         autodiff(Reverse, _apply_instrument!, Duplicated(out, Δout), Duplicated(vis, dvis), Const(J), Const((;)))
#         return NoTangent(), dvis, NoTangent(), NoTangent()
#     end
#     return out, _apply_instrument_pb
# end
