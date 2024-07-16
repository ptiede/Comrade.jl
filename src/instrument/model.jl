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

Base.show(io::IO, mime::MIME"text/plain", m::IdealInstrumentModel) = printstyled(io, "IdealInstrumentModel"; color=:light_cyan, bold=true)

apply_instrument(vis, ::IdealInstrumentModel, x) = vis


struct InstrumentModel{J<:AbstractJonesMatrix, PI, P<:PolBasis} <: AbstractInstrumentModel
    jones::J
    prior::PI
    refbasis::P
end

function Base.show(io::IO, ::MIME"text/plain", m::InstrumentModel)
    printstyled(io, "InstrumentModel"; bold=true, color=:light_cyan)
    println(io)
    T = typeof(m.jones)
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

# Site lookup is const so we add a method so we can signal
# to Enzyme that it is not differentiable.
sitelookup(x::ObservedInstrumentModel) = x.bsitelookup
Enzyme.EnzymeRules.inactive(::typeof(sitelookup), args...) = nothing

function Base.show(io::IO, mime::MIME"text/plain", m::ObservedInstrumentModel)
    printstyled(io, "ObservedInstrumentModel"; bold=true, color=:light_cyan)
    println(io)
    T = typeof(m.instrument)
    ST = split(split(" $T", '{')[1], ".")[end]
    println(io, "  with Jones: ", ST)
    print(io, "  with reference basis: ", m.refbasis)
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
julia> G = SingleStokesGain(x->exp(x.lg + 1im*x.pg))
julia> intprior = (lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1))),
            pg = ArrayPrior(IIDSitePrior(ScanSeg(), DiagVonMises(0.0, inv(π^2))))
            )

julia> intm = InstrumentModel(G, intprior)
```

A standard polarized example is
```julia-repl
julia> G = JonesG() do
        gR = exp.(x.lgr + 1im*x.gpr)
        gL = gr*exp.(x.lgrat + 1im*x.gprat)
        return gR, gL
    end
julia> D = JonesD() do
        dR = complex.(x.dRre, x.dRim)
        dL = complex.(x.dLre, x.dLim)
        return gR, gL
    end
julia> R = JonesR()
julia> J = JonesSandwich(G, D, R)
julia> intprior = (lgr = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)),
                    gpr = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))),
                    lgrat = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)),
                    gprat = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))),
                    dRre = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.1)),
                    dRim = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.1)),
                    dLre = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.1)),
                    dLim = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.1))
                    )
julia> intm = InstrumentModel(J, intprior)
```

which construct the gain matrix from R and ratios, and D is the small leakage matrix. [`JonesR`](@ref)
is the *response matrix* that controls how the site responds to the ideal visibility in the reference
basis.
"""
function InstrumentModel(jones::AbstractJonesMatrix, prior::NamedTuple{N, <:NTuple{M, ArrayPrior}}; refbasis = CirBasis()) where {N, M}
    return InstrumentModel(jones, prior, refbasis)
end

function InstrumentModel(jones::JonesR; refbasis=CirBasis())
    return InstrumentModel(jones, NamedTuple(), refbasis)
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
        i1 = findall(x->(t∈x[1])&&(x[2]==s1), tsf)
        i2 = findall(x->(t∈x[1])&&(x[2]==s2), tsf)
        length(i1) > 1 && throw(AssertionError("Multiple indices found for $t, $((s1)) in SiteArray"))
        length(i2) > 1 && throw(AssertionError("Multiple indices found for $t, $((s2)) in SiteArray"))
        isnothing(i1) && throw(AssertionError("$t, $f, $((s1)) not found in SiteArray"))
        isnothing(i2) && throw(AssertionError("$t, $f, $((s2)) not found in SiteArray"))
        ind1[i] = i1[begin]
        ind2[i] = i2[begin]
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
    _apply_instrument!(baseimage(vout), baseimage(vis), J, x.instrument)
    return vout
end

function apply_instrument(vis, J::ObservedInstrumentModel{<:Union{JonesR, JonesF}}, x)
    vout = intout(vis)
    _apply_instrument!(baseimage(vout), baseimage(vis), J, (;))
    return vout
end


function _apply_instrument!(vout, vis, J::ObservedInstrumentModel, xint)
    # @inbounds for i in eachindex(vout, vis)
    #     vout[i] = apply_jones(vis[i], i, J, xint)
    # end
    vout .= apply_jones.(vis, eachindex(vis), Ref(J), Ref(xint))
    return nothing
end

@inline get_indices(bsitemaps, index, ::Val{1}) = map(x->getindex(x.indices_1, index), bsitemaps)
@inline get_indices(bsitemaps, index, ::Val{2}) = map(x->getindex(x.indices_2, index), bsitemaps)
@inline get_params(x::NamedTuple{N}, indices::NamedTuple{N}) where {N} = NamedTuple{N}(map((xx, ii)->getindex(xx, ii), x, indices))

# We need this because Enzyme seems to crash when generating code for this
# TODO try to find MWE and post to Enzyme.jl
Enzyme.EnzymeRules.inactive(::typeof(get_indices), args...) = nothing

@inline function build_jones(index::Int, J::ObservedInstrumentModel, x, ::Val{N}) where N
    indices = get_indices(sitelookup(J), index, Val(N))
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




function ChainRulesCore.rrule(::typeof(apply_instrument), vis, J::ObservedInstrumentModel, x)
    out = apply_instrument(vis, J, x)
    px = ProjectTo(x)
    function _apply_instrument_pb(Δ)
        bvis = baseimage(vis)
        bout = baseimage(out)
        Δout = similar(bout)
        Δout .= unthunk(Δ)
        xi = x.instrument
        dx = ntzero(xi)
        dvis = zero(bvis)
        autodiff(Reverse, _apply_instrument!, Const, Duplicated(bout, Δout), Duplicated(bvis, dvis), Const(J), Duplicated(xi, dx))
        return NoTangent(), UnstructuredMap(dvis, axisdims(vis)), NoTangent(), px((;instrument = dx))
    end
    return out, _apply_instrument_pb
end

function ChainRulesCore.rrule(::typeof(apply_instrument), vis, J::ObservedInstrumentModel{<:Union{JonesR, JonesF}}, x)
    out = apply_instrument(vis, J, x)
    function _apply_instrument_pb(Δ)
        Δout = similar(out)
        Δout .= unthunk(Δ)
        dvis = zero(vis)
        autodiff(Reverse, _apply_instrument!, Duplicated(out, Δout), Duplicated(vis, dvis), Const(J), Const((;)))
        return NoTangent(), dvis, NoTangent(), NoTangent()
    end
    return out, _apply_instrument_pb
end
