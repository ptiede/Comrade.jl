export SingleStokesGain, JonesG, JonesD, JonesF, JonesR, GenericJones,
       JonesModel

abstract type AbstractJonesMatrix end
@inline jonesmatrix(mat::AbstractJonesMatrix, params, visindex, site) = construct_jones(mat, param_map(mat, params), visindex, site)
@inline param_map(mat::AbstractJonesMatrix, x) = mat.param_map(x)
preallocate_jones(g::AbstractJonesMatrix, array, refbasis) = g

struct SingleStokesGain{F} <: AbstractJonesMatrix
    param_map::F
end
@inline contruct_jones(::SingleStokesGain, x, index, site) = x

struct JonesG{F} <: AbstractJonesMatrix
    param_map::F
end
construct_jones(::JonesG, x::NTuple{2, T}, index, site) where {T} = Diagonal(SVector{2, T}(x))


struct JonesD{F} <: AbstractJonesMatrix
    param_map::F
end
construct_jones(::JonesD, x::NTuple{2, T}, index, site) where {T} = SMatrix{2, 2, T, 4}(1, x[1], x[2], 1)


"""
    GenericJones

Construct a generic dense jones matrix with four parameterized elements
"""
struct GenericJones{F} <: AbstractJonesMatrix
    param_map::F
end
construct_jones(::GenericJones, x::NTuple{4, T}, index, site) where {T} = SMatrix{2, 2, T, 4}(x[1], x[2], x[3], x[4])

struct JonesF{M} <: AbstractJonesMatrix
    matrices::M
end
JonesF() = JonesF(nothing)
construct_jones(J::JonesF, x, index, ::Val{M}) where {M} = J.matrices[index][M]
function preallocate_jones(::JonesF, array::ArrayConfiguration)
    field_rotations = build_frs(array)
    return field_rotations
end

Base.@kwdef struct JonesR{M} <: AbstractJonesMatrix
    matrices::M = nothing
    add_fr::Bool = true
end
construct_jones(J::JonesR, x, index, ::Val{M}) where {M} = J.matrices[index][M]
function preallocate_jones(J::JonesR, array::ArrayConfiguration, ref)
    T1 = StructArray(map(x -> basis_transform(ref, x[1]), array.data.polbasis))
    T2 = StructArray(map(x -> basis_transform(ref, x[2]), array.data.polbasis))
    Tcirc1 = StructArray(map(x -> basis_transform(CirBasis(), x[1]), array.data.polbasis))
    Tcirc2 = StructArray(map(x -> basis_transform(CirBasis(), x[2]), array.data.polbasis))
    if J.add_fr
        field_rotations = build_frs(array)
        @. T1 .= Tcirc1*field_rotations.m1*adjoint(Tcirc1)*T1
        @. T2 .= Tcirc2*field_rotations.m2*adjoint(Tcirc2)*T2
    end
    return JonesPairs(T1, T2)

end


struct JonesModel{J, M} <: AbstractJonesMatrix
    jones_map::J
    matrices::M
end

function JonesModel(map, matrices::AbstractJonesMatrix...)
    return JonesModel(map, matrices)
end

param_map(j::JonesModel, x) = map(j->param_map(j, x), j.matrices)
function preallocate_jones(J::JonesModel, array::ArrayConfiguration)
    m2 = map(x->preallocate_jones(x, array), J.matrices, J.refbasis)
    return JonesModel(J.jones_map, m2, J.refbasis)
end
