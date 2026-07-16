export SingleStokesGain, JonesG, JonesD, JonesF, JonesR, GenericJones,
    JonesSandwich, JonesPoint, forward_jones

abstract type AbstractJonesMatrix end

"""
    JonesPoint

The metadata of one gain application point, passed as the *second* argument of every
`param_map`. Fields:

  - `Ti`: the datum's time in UTC hours
  - `Fr`: the datum's frequency in Hz. This is always the frequency of the *datum* being
    corrupted (never a channel center), so frequency-dependent gain laws compose with
    gains shared across channels (`freqseg = FullBand()`).
  - `site`: the antenna's site `Symbol`. Under Reactant compilation this is `nothing`
    (symbols cannot be gathered inside a compiled loop) — encode site dependence in the
    priors, not in `param_map` branches, if you intend to compile.
  - `vind`: the linear index of the datum in the visibility table.

The antenna slot (1 or 2 of the baseline) is the compile-time parameter `N`, retrievable
with `antslot(p)`.

# Example

```julia
G = SingleStokesGain() do x, p
    lg = x.lg0 + x.α * log(p.Fr / 230.0e9)
    return exp(lg + 1im * x.gp)
end
```
"""
struct JonesPoint{N, T, F, S, I}
    Ti::T
    Fr::F
    site::S
    vind::I
end

@inline JonesPoint{N}(Ti::T, Fr::F, site::S, vind::I) where {N, T, F, S, I} =
    JonesPoint{N, T, F, S, I}(Ti, Fr, site, vind)

"""
    antslot(p::JonesPoint)

The antenna slot of the gain point: `1` for the first site of the baseline, `2` for the
second. A compile-time constant.
"""
antslot(::JonesPoint{N}) where {N} = N

@inline jonesmatrix(mat::AbstractJonesMatrix, params, p::JonesPoint) =
    construct_jones(mat, param_map(mat, params, p), p)
@inline param_map(mat::AbstractJonesMatrix, x, p) = mat.param_map(x, p)
preallocate_jones(g::AbstractJonesMatrix, array, refbasis) = g

"""
    SingleStokesGain(param_map)

Construct a gain term that is applicable to a single measured visibility. This is
useful for pure stokes I modeling. The `param_map` is a function `(x, p) -> gain` where
`x` is a named tuple of the sampled parameter values at the gain point and `p` is the
point's [`JonesPoint`](@ref) metadata (time, frequency, site).
The return value of the `param_map` should be a single number or complex gain.

## Example
```julia
G = SingleStokesGain((x, p) -> exp(x.lg + 1im * x.gp))
```

"""
struct SingleStokesGain{F} <: AbstractJonesMatrix
    param_map::F
end
@inline construct_jones(::SingleStokesGain, x, ::JonesPoint) = x

"""
    JonesG(param_map)

Describes a gain Jones matrix with layout

  g1 0
  0  g2

where `g1` and `g2` are the gains for first and second feed of the telescope.

The `param_map` is a function `(x, p) -> (g1, g2)` where `x` is a named tuple of the
sampled parameter values at the gain point and `p` is the point's [`JonesPoint`](@ref)
metadata (time, frequency, site).
The return value of the `param_map` should be a two element tuple where the first element
is the complex gain `g1` and the second element is the complex gain `g2`.

## Example
```julia
G = JonesG() do x, p
    g1 = exp(complex(x.lg1, x.gp1))
    g2 = g1*exp(complex(x.lgratio, x.gpratio))
    return g1, g2
end
```
"""
struct JonesG{F} <: AbstractJonesMatrix
    param_map::F
end
@inline construct_jones(::JonesG, x::NTuple{2, T}, ::JonesPoint) where {T} = SMatrix{2, 2, T, 4}(x[1], zero(T), zero(T), x[2])


"""
    JonesD(param_map)

Describes a leakage Jones matrix with layout

  1 d1
  d2 1

where `d1` and `d2` are the d-terms for first and second feed of the telescope.

The `param_map` is a function `(x, p) -> (d1, d2)` where `x` is a named tuple of the
sampled parameter values at the gain point and `p` is the point's [`JonesPoint`](@ref)
metadata (time, frequency, site).
The return value of the `param_map` should be a two element tuple where the first element
is the complex d-term `d1` and the second element is the complex d-term `d2`.

## Example
```julia
D = JonesD() do x, p
    d1 = complex(x.d1real, x.d1imag)
    d2 = complex(x.d2real, x.d2imag)
    return d1, d2
end
```
"""
struct JonesD{F} <: AbstractJonesMatrix
    param_map::F
end
Base.@propagate_inbounds construct_jones(::JonesD, x::NTuple{2, T}, ::JonesPoint) where {T} = SMatrix{2, 2, T, 4}(1, x[2], x[1], 1)


"""
    GenericJones(param_map)

Construct a generic dense jones matrix with four parameterized elements.

The `param_map` is a function `(x, p) -> (j11, j21, j12, j22)` where `x` is a named tuple
of the sampled parameter values at the gain point and `p` is the point's
[`JonesPoint`](@ref) metadata (time, frequency, site).
The return value of the `param_map` should be a four element tuple where the elements are
the entries of the jones matrix in column major order.

## Example
```julia
J = GenericJones() do x, p
    return x.j11, x.j21, x.j12, x.j22
end
```

"""
struct GenericJones{F} <: AbstractJonesMatrix
    param_map::F
end
@inline construct_jones(::GenericJones, x::NTuple{4, T}, ::JonesPoint) where {T} = SMatrix{2, 2, T, 4}(x[1], x[2], x[3], x[4])


struct JonesConst{M} <: AbstractJonesMatrix
    matrices::M
end
JonesConst(m1, m2) = JonesConst((m1, m2))
@inline construct_jones(J::JonesConst, x, p::JonesPoint{N}) where {N} = @inbounds(rgetindex(J.matrices[N], p.vind))
param_map(::JonesConst, x, p) = x


"""
    JonesF(;add_fr=true)

The **feed rotation** Jones matrix. This matrix describes the orientation of the feeds of the
telescope.

This Jones matrix has no parameters so it doesn't accept
a `param_map`. The `add_fr` argument is a boolean that specifies if feed rotation should
be included.
"""
struct JonesF{M} <: AbstractJonesMatrix
    matrices::M
end
JonesF() = JonesF(nothing)
@inline construct_jones(J::JonesF, x, p::JonesPoint{N}) where {N} = @inbounds J.matrices[p.vind][N]
param_map(::JonesF, x, p) = x
function preallocate_jones(::JonesF, array::AbstractArrayConfiguration, ref)
    field_rotations = build_feedrotation(array)
    return JonesConst(field_rotations[1], field_rotations[2])
end


"""
    JonesR(;add_fr=true)

The **response** Jones matrix. This is the reponse the telescope has to the incoming
electric field, if the telescope was ideal. If `add_fr=true` then feed rotation are included.

This Jones matrix has no parameters so it doesn't accept
a `param_map`. The `add_fr` argument is a boolean that specifies if feed rotation should
be included.
"""
Base.@kwdef struct JonesR{M} <: AbstractJonesMatrix
    matrices::M = nothing
    add_fr::Bool = true
end
@inline construct_jones(J::JonesR, x, p::JonesPoint{N}) where {N} = @inbounds rgetindex(J.matrices[N], p.vind)
param_map(::JonesR, x, p) = x

function preallocate_jones(J::JonesR, array::AbstractArrayConfiguration, ref)
    T1 = StructArray(map(x -> basis_transform(ref, x[1]), array[:polbasis]))
    T2 = StructArray(map(x -> basis_transform(ref, x[2]), array[:polbasis]))
    Tcirc1 = StructArray(map(x -> basis_transform(CirBasis(), x[1]), array[:polbasis]))
    Tcirc2 = StructArray(map(x -> basis_transform(CirBasis(), x[2]), array[:polbasis]))
    if J.add_fr
        field_rotations = build_feedrotation(array)
        @. T1 .= Tcirc1 * field_rotations[1] * adjoint(Tcirc1) * T1
        @. T2 .= Tcirc2 * field_rotations[2] * adjoint(Tcirc2) * T2
    end
    return JonesConst(T1, T2)

end


struct JonesSandwich{J, M} <: AbstractJonesMatrix
    jones_map::J
    matrices::M
end

"""
    JonesSandwich([decomp_function=*,] matrices::AbstractJonesMatrix...)

Constructs a Jones matrix that is the results combining multiple Jones matrices together.
The specific composition is determined by the `decomp_function`. For example if the
decomp function is `*` then the matrices are multiplied together, if it is `+` then they
are added.


## Examples
```julia
G = JonesG((x, p)->(x.gR, x.gL)) # Gain matrix
D = JonesD((x, p)->(x.dR, x.dL)) # leakage matrix
F = JonesF()                # Feed rotation matrix

J = JonesSandwich(*, G, D, F) # Construct the full Jones matrix as G*D*F

# Or if you want to include FR calibration
J = JonesSandwich(G, D, F) do g, d, f
    return adjoint(f)*g*d*f
end
```
"""
function JonesSandwich(map, matrices::AbstractJonesMatrix...)
    return JonesSandwich(splat(map), matrices)
end

function JonesSandwich(matrices::AbstractJonesMatrix...)
    return JonesSandwich(*, matrices...)
end

@inline function jonesmatrix(J::JonesSandwich, x, p::JonesPoint)
    return J.jones_map(map(m -> construct_jones(m, param_map(m, x, p), p), J.matrices))
end

function preallocate_jones(J::JonesSandwich, array::AbstractArrayConfiguration, refbasis = CirBasis())
    m2 = map(x -> preallocate_jones(x, array, refbasis), J.matrices)
    return JonesSandwich(J.jones_map, m2)
end

"""
    forward_jones(J::AbstractJonesMatrix, xs::NamedTuple)

Construct the Jones matrices of the gain model `J` at every point of the shared grid of
the parameters `xs` (a NamedTuple of `SiteArray`s, one per parameter). Every element of
`xs` must be on the *same* (time, frequency, site) grid — the evaluation is pointwise,
with no matching involved. For parameters with mixed segmentations use
[`instrumentmodel(post, θ)`](@ref), which evaluates through the observed instrument
model's exact per-datum index maps instead.
"""
function forward_jones(v::AbstractJonesMatrix, xs::NamedTuple{N}) where {N}
    x1 = first(values(xs))
    for x in values(xs)
        (times(x) == times(x1) && frequencies(x) == frequencies(x1) && sites(x) == sites(x1)) ||
            throw(
            ArgumentError(
                "forward_jones requires every element of `xs` to share one " *
                    "(time, frequency, site) grid. Parameters with mixed segmentations " *
                    "should be evaluated through the observed instrument model, e.g. " *
                    "`instrumentmodel(post, θ)`."
            )
        )
    end
    T = times(x1)
    F = frequencies(x1)
    S = sites(x1)
    vs = [
        jonesmatrix(
                v, map(Base.Fix2(getindex, i), xs),
                JonesPoint{1}(float(T[i]), float(F[i]), S[i], i)
            )
            for i in eachindex(parent(x1))
    ]
    return SiteArray(vs, T, F, S)
end

# the mixed-segmentation form, `forward_jones(::ObservedInstrumentModel, xs)`, lives in
# model.jl with the observed model and its per-datum index maps
