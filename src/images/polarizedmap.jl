export PolarizedMap, stokes_parameter, StokesVector, CoherencyMatrix, evpa, m̆

struct StokesVector{T} <: FieldVector{4,T}
    i::T
    q::T
    u::T
    v::T
end

struct CoherencyMatrix{T} <: FieldMatrix{2,2,T}
    rr::T
    lr::T
    rl::T
    ll::T
end


@inline function Base.convert(::Type{CoherencyMatrix}, p::StokesVector)
    rr = p.i + p.v
    ll = p.i - p.v
    rl = p.q + 1im*p.u
    lr = p.q - 1im*p.u
    return CoherencyMatrix(rr, lr, rl, ll)
end

@inline function Base.convert(::Type{StokesVector}, p::CoherencyMatrix)
    i = (p.rr + p.ll)/2
    v = (p.rr - p.ll)/2
    q = (p.rl + p.lr)/2
    u = (p.rl - p.lr)/(2im)
    return CircularMatrix(i, q, u, v)
end


m̆(m::StokesVector) = m.q + 1im*m.u
m̆(m::CoherencyMatrix) = m.rl


evpa(m::StokesVector) = 1/2*atan2(m.U, m.Q)
evpa(m::CoherencyMatrix) = evpa(convert(StokesVector, m))



struct PolarizedMap{SI<:AbstractIntensityMap,
                    SQ<:AbstractIntensityMap,
                    SU<:AbstractIntensityMap,
                    SV<:AbstractIntensityMap} <: AbstractPolarizedMap{SI,SQ,SU,SV}
    I::SI
    Q::SQ
    U::SU
    V::SV
    function PolarizedMap(I::SI,Q::SQ,U::SU,V::SV) where {SI, SQ, SU, SV}
        @assert size(I) == size(Q) == size(U) == size(V) "Image sizes must be equal in polarized map"
        @assert fov(I) == fov(Q) == fov(U) == fov(V) "Image fov must be equal in polarized map"
        new{SI,SQ,SU,SV}(I,Q,U,V)
    end
end


Base.Base.@propagate_inbounds function Base.getindex(pimg::PolarizedMap, i...)
return StokesVector(pimg.I[i...], pimg.Q[i...], pimg.U[i...], pimg.V[i...])
end

@inline stokes_parameter(pimg::PolarizedMap, p::Symbol) = getproperty(pimg, p)
