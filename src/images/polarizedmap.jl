export PolarizedMap, stokes_parameter

struct PolarizedMap{SI<:AbstractIntensityMap,
    SQ<:AbstractIntensityMap,
    SU<:AbstractIntensityMap,
    SV<:AbstractIntensityMap} <: AbstractPolarizedMap
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
