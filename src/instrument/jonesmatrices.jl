export Jones, PerStation, gains, leakage, rotation, forward_jones

"""
    AbstractJonesMatrix

An instrument term in the local ideal feed basis of a station: per station cell, a 2×2 matrix,
or a number for a term that acts on both feeds alike. Each term owns the priors of the
parameters it reads. Build one with [`Jones`](@ref) and combine terms with `*`.
"""
abstract type AbstractJonesMatrix end

"""
    Jones(f; name = prior, ...)

An instrument term with the parameters `name`, each with its prior (an [`ArrayPrior`](@ref)).
`f(x, c)` takes the instrument parameters `x` of one station cell and the cell's coordinates
`c` (`c.ant`, the station's index in the station table; `c.Ti`, hours since the epoch; `c.Fr`,
Hz; `c.elevation`, `c.parallactic`, `c.azimuth`, radians; `c.ra`, `c.dec`, the J2000 direction
of the observed source, radians) and returns the term's 2×2 matrix in
the station's local feed basis, or a number, a gain on both feeds (the form for Stokes I
data). [`gains`](@ref), [`leakage`](@ref) and [`rotation`](@ref) build the usual matrices; any
other optic is the matrix written out. A parameter exists only at the stations whose term
reads it (see [`PerStation`](@ref)), and two terms that read the same parameter must give it
the same prior.

```julia
G = Jones(
    (x, c) -> gains(exp(complex(x.lgR, x.gpR)), exp(complex(x.lgL, x.gpL)));
    lgR = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1))),
    gpR = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
    lgL = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1))),
    gpL = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
)
# a quarter-wave plate whose retardance and angle are fit
Q = Jones(; θ = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.1))),
            δ = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.1)))) do x, c
    R = rotation(x.θ, LinBasis())
    return R' * SMatrix{2, 2}(1, 0, 0, cis(π / 2 + x.δ)) * R
end
J = G * Q
```
"""
struct Jones{F, P <: NamedTuple} <: AbstractJonesMatrix
    param_map::F
    priors::P
end
Jones(f; priors...) = Jones(f, NamedTuple(priors))
@inline jonesmatrix(J::Jones, x, c) = J.param_map(x, c)

# The priors of the parameters a term reads, by name.
_priors(J::Jones) = J.priors
_priors(::AbstractJonesMatrix) = (;)

"""
    gains(g1, g2)

The gain matrix `[g1 0; 0 g2]`: the complex gains of the first and second feed.
"""
@inline function gains(g1, g2)
    z = zero(promote_type(typeof(g1), typeof(g2)))
    return SMatrix{2, 2}(g1, z, z, g2)
end

"""
    leakage(d1, d2)

The leakage matrix `[1 d1; d2 1]`: `d1` is the leakage of the second feed's polarization into
the first feed, `d2` of the first into the second.
"""
@inline function leakage(d1, d2)
    o = one(promote_type(typeof(d1), typeof(d2)))
    return SMatrix{2, 2}(o, d2, d1, o)
end

"""
    rotation(φ, basis::PolBasis)

The rotation of the electric field by the angle `φ` (radians) against the sky, as a 2×2 matrix
in the feed basis `basis`: `diag(e^{-iφ}, e^{iφ})` in `CirBasis()`, and that matrix transformed
to any other basis. A station's feed rotation (see [`feed_rotation`](@ref)) acts on its feeds
this way.
"""
@inline function rotation(φ, basis::PolBasis)
    Tc = basis_transform(CirBasis(), basis)
    return Tc * SMatrix{2, 2}(cis(-φ), zero(cis(φ)), zero(cis(φ)), cis(φ)) * Tc'
end

"""
    J1 * J2

The product of the instrument terms `J1` and `J2`, evaluated per station cell as the matrix
product in that order, e.g. `G * D` for gains times leakage.
"""
struct JonesProduct{M <: Tuple} <: AbstractJonesMatrix
    matrices::M
end

_factors(J::JonesProduct) = J.matrices
_factors(J::AbstractJonesMatrix) = (J,)
Base.:*(a::AbstractJonesMatrix, b::AbstractJonesMatrix) = JonesProduct((_factors(a)..., _factors(b)...))

@inline jonesmatrix(J::JonesProduct, x, c) = _mulall(map(m -> jonesmatrix(m, x, c), J.matrices)...)

# The product unrolled at compile time: `prod` over a tuple is a loop that Enzyme reverses with
# runtime dispatch.
@inline _mulall(a) = a
@inline _mulall(a, b, rest...) = _mulall(a * b, rest...)

"""
    PerStation(default; station = term, ...)

Instrument terms assigned per station: each named station uses its `term`, every other station
`default`, e.g. `PerStation(G * D; AA = G * D * Q)` for a quarter-wave plate at AA only. Every
term must give values of one type (all 2×2 matrices, or all numbers). A parameter exists only
at the stations whose term reads it: above, the parameters of `Q` exist only at AA, and a
station whose term leaves out `D` has no leakage parameters.
"""
struct PerStation{D <: AbstractJonesMatrix, O <: NamedTuple} <: AbstractJonesMatrix
    default::D
    overrides::O
end

function PerStation(default::AbstractJonesMatrix; kwargs...)
    o = NamedTuple(kwargs)
    all(t -> t isa AbstractJonesMatrix, values(o)) ||
        throw(ArgumentError("PerStation takes an instrument term for each station"))
    any(t -> t isa PerStation, values(o)) && throw(ArgumentError("PerStation terms cannot be nested"))
    default isa PerStation && throw(ArgumentError("PerStation terms cannot be nested"))
    return PerStation(default, o)
end

Base.:*(::PerStation, ::AbstractJonesMatrix) = throw(ArgumentError("a PerStation model is the whole instrument; multiply the terms inside it"))
Base.:*(::AbstractJonesMatrix, ::PerStation) = throw(ArgumentError("a PerStation model is the whole instrument; multiply the terms inside it"))
Base.:*(::PerStation, ::PerStation) = throw(ArgumentError("a PerStation model is the whole instrument; multiply the terms inside it"))

# The term of station `s`.
stationterm(J::PerStation, s::Symbol) = get(J.overrides, s, J.default)

_priors(J::JonesProduct) = _mergepriors(map(_priors, J.matrices)...)
_priors(J::PerStation) = _mergepriors(_priors(J.default), map(_priors, values(J.overrides))...)

_mergepriors() = (;)
_mergepriors(a::NamedTuple) = a
function _mergepriors(a::NamedTuple, b::NamedTuple, rest::NamedTuple...)
    for k in keys(b)
        haskey(a, k) && getproperty(a, k) !== getproperty(b, k) && throw(
            ArgumentError("two instrument terms give the parameter $k different priors; a parameter shared by terms has one prior")
        )
    end
    return _mergepriors(merge(a, b), rest...)
end

# Whether the term `J` reads the parameter `p`.
_reads(J::AbstractJonesMatrix, p) = haskey(_priors(J), p)

# `J` as the observed model evaluates it: the priors are needed only to build the model's prior.
_withoutpriors(J::Jones) = Jones(J.param_map, (;))
_withoutpriors(J::JonesProduct) = JonesProduct(map(_withoutpriors, J.matrices))
_withoutpriors(J::PerStation) = PerStation(_withoutpriors(J.default), map(_withoutpriors, J.overrides))
_withoutpriors(J::AbstractJonesMatrix) = J
