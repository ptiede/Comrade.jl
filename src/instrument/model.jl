export InstrumentModel

abstract type AbstractInstrumentModel end

"""
    IdealInstrument(array::AbstractArrayConfiguration)

Constructs an ideal instrument that has no corruptions including feed rotations.
"""
struct IdealInstrumentModel <: AbstractInstrumentModel end

Base.show(io::IO, mime::MIME"text/plain", m::IdealInstrumentModel) = printstyled(io, "IdealInstrumentModel"; color=:light_cyan, bold=true)

struct InstrumentModel{J<:AbstractJonesMatrix, PI, P<:PolBasis} <: AbstractInstrumentModel
    jones::J
    prior::PI
    refbasis::P
end

function Base.show(io::IO, ::MIME"text/plain", m::InstrumentModel)
    printstyled(io, "InstrumentModel"; bold=true, color=:light_cyan)
    println(io)
    T = typeof(m.instrument)
    ST = split(split(" $T", '{')[1], ".")[end]
    println(io, "  with Jones: ", ST)
    print(io, "  with reference basis: ", m.refbasis)
end




struct ObservedInstrumentModel{I<:AbstractJonesMatrix, PB<:PolBasis, B} <: AbstractInstrumentModel
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
end

function Base.show(io::IO, mime::MIME"text/plain", m::ObservedInstrumentModel)
    printstyled(io, "ObservedInstrumentModel"; bold=true, color=:light_cyan)
    println(io)
    T = typeof(m.instrument)
    ST = split(split(" $T", '{')[1], ".")[end]
    println(io, "  with Jones: ", ST)
    print(io, "  with reference basis: ", m.refbasis)
end


"""
    InstrumentModel(jones, prior, array; refbasis = CirBasis())

Builds an instrument model using the jones matrix `jones`, with priors `prior` and the
array configuration `array`. The reference basis is `refbasis` and is used to define what
the ideal basis is. Namely, the basis that you have the ideal visibilties to be represented in.
"""
function InstrumentModel(jones::AbstractJonesMatrix, prior::NamedTuple{N, <:NTuple{M, ArrayPrior}}; refbasis = CirBasis()) where {N, M}
    return InstrumentModel(jones, prior, refbasis)
end

function set_array(int::InstrumentModel, array::AbstractArrayConfiguration)
    (;jones, prior, refbasis) = int
    # 1. preallocate and jones matrices
    Jpre = preallocate_jones(jones, array, refbasis)
    # 2. construct the prior with the array you have
    prior_obs = NamedDist(map(x->ObservedArrayPrior(x, array), prior))
    # 3. construct the baseline site map for each prior
    x = rand(prior_obs)
    bsitemaps = map(x->_construct_baselinemap(array, x), x)
    intobs = ObservedInstrumentModel(Jpre, refbasis, bsitemaps)
    return intobs, prior_obs
end

function set_array(int::IdealInstrumentModel, ::AbstractArrayConfiguration)
    return (int, ())
end

struct BaselineSiteLookup{V<:AbstractArray{<:Integer}}
    indices_1::V
    indices_2::V
end

function _construct_baselinemap(array::EHTArrayConfiguration, x::SiteArray)
    T = array[:Ti]
    F = array[:Fr]
    bl = array[:sites]

    return _construct_baselinemap(T, F, bl, x)
end

function _construct_baselinemap(T, F, bl, x::SiteArray)
    tcal = times(x)
    scal = sites(x)
    fcal = frequencies(x)
    tsf = StructArray((tcal, scal, fcal))
    ind1 = similar(T, Int)
    ind2 = similar(T, Int)
    for i in eachindex(T, F, bl, ind1, ind2)
        t = T[i]
        f = F[i]
        s1, s2 = bl[i]
        i1 = findfirst(x->(t∈x[1])&&(x[2]==s1), tsf)
        i2 = findfirst(x->(t∈x[1])&&(x[2]==s2), tsf)
        isnothing(i1) && throw(AssertionError("$t, $f, $((s1)) not found in SiteArray"))
        isnothing(i2) && throw(AssertionError("$t, $f, $((s2)) not found in SiteArray"))
        ind1[i] = i1
        ind2[i] = i2
    end
    BaselineSiteLookup(ind1, ind2)
end


intout(vis::AbstractArray{<:StokesParams{T}}) where {T<:Real} = similar(vis, SMatrix{2,2, Complex{T}, 4})
intout(vis::AbstractArray{T}) where {T<:Real} = similar(vis, Complex{T})
intout(vis::AbstractArray{<:CoherencyMatrix{A,B,T}}) where {A,B,T<:Real} = similar(vis, SMatrix{2,2, Complex{T}, 4})

intout(vis::AbstractArray{<:StokesParams{T}}) where {T<:Complex} = similar(vis, SMatrix{2,2, T, 4})
intout(vis::AbstractArray{T}) where {T<:Complex} = similar(vis, T)
intout(vis::AbstractArray{<:CoherencyMatrix{A,B,T}}) where {A,B,T<:Complex} = similar(vis, SMatrix{2,2, T, 4})


function apply_instrument(vis, J::ObservedInstrumentModel, x)
    vout = intout(vis)
    _apply_instrument!(vout, vis, J, x.instrument)
    return vout
end

function _apply_instrument!(vout, vis, J::ObservedInstrumentModel, xint)
    @inbounds for i in eachindex(vout, vis)
        vout[i] = apply_jones(vis[i], i, J, xint)
    end
    # vout .= apply_jones.(vis, eachindex(vis), Ref(J), Ref(x))
    return nothing
end

@inline get_indices(bsitemaps, index, ::Val{1}) = map(x->getindex(x.indices_1, index), bsitemaps)
@inline get_indices(bsitemaps, index, ::Val{2}) = map(x->getindex(x.indices_2, index), bsitemaps)
@inline get_params(x::NamedTuple{N}, indices::NamedTuple{N}) where {N} = NamedTuple{N}(map((xx, ii)->getindex(xx, ii), x, indices))
@inline function build_jones(index::Int, J::ObservedInstrumentModel, x, ::Val{N}) where N
    indices = get_indices(J.bsitelookup, index, Val(N))
    params = get_params(x, indices)
    return jonesmatrix(J.instrument, params, index, Val(N))
end


@inline function apply_jones(v, index::Int, J::ObservedInstrumentModel, x)
    j1 = build_jones(index, J, x, Val(1))
    j2 = build_jones(index, J, x, Val(2))
    vout =  _apply_jones(v, j1, j2, J.refbasis)
    return vout
end


@inline _apply_jones(v::Number, j1, j2, ::B) where {B} = j1*v*conj(j2)
@inline _apply_jones(v::CoherencyMatrix, j1, j2, ::B) where {B} = j1*CoherencyMatrix{B,B}(v)*j2'
@inline _apply_jones(v::StokesParams, j1, j2, ::B) where {B} = j1*CoherencyMatrix{B,B}(v)*j2'



apply_instrument(vis, ::IdealInstrumentModel, x) = vis

function ChainRulesCore.rrule(::typeof(apply_instrument), vis, J::ObservedInstrumentModel, x)
    out = apply_instrument(vis, J, x)
    function _apply_instrument_pb(Δ)
        Δout = similar(out)
        Δout .= unthunk(Δ)
        xi = x.instrument
        dx = ntzero(xi)
        dvis = zero(vis)
        autodiff(Reverse, _apply_instrument!, Duplicated(out, Δout), Duplicated(vis, dvis), Const(J), Duplicated(xi, dx))
        return NoTangent(), dvis, NoTangent(), Tangent{typeof(x)}(;sky = ZeroTangent(), instrument = dx)
    end
    return out, _apply_instrument_pb
end
