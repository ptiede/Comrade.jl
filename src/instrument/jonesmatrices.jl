export SingleStokesGain, JonesG, JonesD, JonesF, JonesR, GenericJones,
       JonesSandwich

abstract type AbstractJonesMatrix end
@inline jonesmatrix(mat::AbstractJonesMatrix, params, visindex, site) = construct_jones(mat, param_map(mat, params), visindex, site)
@inline param_map(mat::AbstractJonesMatrix, x) = mat.param_map(x)
preallocate_jones(g::AbstractJonesMatrix, array, refbasis) = g

"""
    SingleStokesGain(param_map)

Construct a gain term that is applicable to a single measured visibility. This is
useful for pure stokes I modeling. The `param_map` is a function that maps the
set of parameters to the gain term.

The arguments of `param_map` is typically a named tuple, where some of the elements
are the parameter values needed to compute the gain.
The return value of the `param_map` should be a single number or complex gain.

## Example
```julia
G = SingleStokesGain(x->exp(x.lg + xp.gp))
end
"""
struct SingleStokesGain{F} <: AbstractJonesMatrix
    param_map::F
end
construct_jones(::SingleStokesGain, x, index, site) = x

"""
    JonesG(param_map)

Describes a gain Jones matrix with layout

  g1 0
  0  g2

where `g1` and `g2` are the gains for first and second feed of the telescope.

The `param_map` is a function that maps the set of parameters to the gain term.
The arguments of `param_map` is typically a named tuple, where some of the elements
are the parameter values needed to compute the gain.
The return value of the `param_map` should be a two element tuple where the first element
is the complex gain `g1` and the second element is the complex gain `g2`.

## Example
```julia
G = JonesG() do
    g1 = exp(x.lg1 + 1im.*x.gp1)
    g2 = g1*exp(x.lgratio + 1im.*x.gpratio)
    return g1, g2
end
"""
struct JonesG{F} <: AbstractJonesMatrix
    param_map::F
end
construct_jones(::JonesG, x::NTuple{2, T}, index, site) where {T} = Diagonal(SVector{2, T}(x))


"""
    JonesG(param_map)

Describes a leakage Jones matrix with layout

  1 d1
  d2 1

where `d1` and `d2` are the d-terms for first and second feed of the telescope.

The `param_map` is a function that maps the set of parameters to the gain term.
The arguments of `param_map` is typically a named tuple, where some of the elements
are the parameter values needed to compute the gain.
The return value of the `param_map` should be a two element tuple where the first element
is the complex d-term `d1` and the second element is the complex d-term `d2`.

## Example
```julia
D = JonesD() do
    d1 = complex(x.d1real, x.d1imag)
    d2 = complex(x.d2real, x.d2imag)
    return d1, d2
end
"""
struct JonesD{F} <: AbstractJonesMatrix
    param_map::F
end
construct_jones(::JonesD, x::NTuple{2, T}, index, site) where {T} = SMatrix{2, 2, T, 4}(1, x[2], x[1], 1)


"""
    GenericJones(param_map)

Construct a generic dense jones matrix with four parameterized elements.

The `param_map` is a function that maps the set of parameters to the gain term.
The arguments of `param_map` is typically a named tuple, where some of the elements
are the parameter values needed to compute the gain.
The return value of the `param_map` should be a four element tuple where the elements are
the entries of the jones matrix in column major order.

## Example
```julia
J = GenericJones() do
    return x.j11, x.j21, x.j12, x.j22
end

"""
struct GenericJones{F} <: AbstractJonesMatrix
    param_map::F
end
construct_jones(::GenericJones, x::NTuple{4, T}, index, site) where {T} = SMatrix{2, 2, T, 4}(x[1], x[2], x[3], x[4])

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
construct_jones(J::JonesF, x, index, ::Val{M}) where {M} = J.matrices[index][M]
param_map(::JonesF, x) = x
function preallocate_jones(::JonesF, array::AbstractArrayConfiguration)
    field_rotations = build_frs(array)
    return JonesF(field_rotations)
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
construct_jones(J::JonesR, x, index, ::Val{M}) where {M} = J.matrices[M][index]
param_map(::JonesR, x) = x

function preallocate_jones(J::JonesR, array::AbstractArrayConfiguration, ref)
    T1 = StructArray(map(x -> basis_transform(ref, x[1]), array[:polbasis]))
    T2 = StructArray(map(x -> basis_transform(ref, x[2]), array[:polbasis]))
    Tcirc1 = StructArray(map(x -> basis_transform(CirBasis(), x[1]), array[:polbasis]))
    Tcirc2 = StructArray(map(x -> basis_transform(CirBasis(), x[2]), array[:polbasis]))
    if J.add_fr
        field_rotations = build_feedrotation(array)
        @. T1 .= Tcirc1*field_rotations[1]*adjoint(Tcirc1)*T1
        @. T2 .= Tcirc2*field_rotations[2]*adjoint(Tcirc2)*T2
    end
    return JonesR((T1, T2), J.add_fr)

end


struct JonesSandwich{J, M} <: AbstractJonesMatrix
    jones_map::J
    matrices::M
end

"""
    JonesSandwich([decomp_function=splat(*),] matrices::AbstractJonesMatrix...)

Constructs a Jones matrix that is the results combining multiple Jones matrices together.
The specific composition is determined by the `decomp_function`. For example if the
decomp function is `*` then the matrices are multiplied together, if it is `+` then they
are added.


## Examples
```julia
G = JonesG(x->(x.gR, x.gL)) # Gain matrix
D = JonesD(x->(x.dR, x.dL)) # leakage matrix
F = JonesF()                # Feed rotation matrix

J = JonesSandwich(splat(*), G, D, F) # Construct the full Jones matrix as G*D*F

# Or if you want to include FR calibration
J = JonesSandwich(G, D, F) do g, d, f
    return adjoint(f)g*d*f
end
```
"""
function JonesSandwich(map, matrices::AbstractJonesMatrix...)
    return JonesSandwich(map, matrices)
end

function JonesSandwich(matrices::AbstractJonesMatrix...)
    return JonesSandwich(splat(*), matrices...)
end

function jonesmatrix(J::JonesSandwich, x, index, site)
    return J.jones_map(map(m->construct_jones(m, param_map(m, x), index, site), J.matrices))
end

function preallocate_jones(J::JonesSandwich, array::AbstractArrayConfiguration, refbasis=CirBasis())
    m2 = map(x->preallocate_jones(x, array, refbasis), J.matrices)
    return JonesSandwich(J.jones_map, m2)
end
