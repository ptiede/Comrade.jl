"""
    $(TYPEDEF)
Adds two models together to create composite models. Note that
I may change this in the future so make it easier on the compiler,
i.e. make the composite model a fancy model vector and heap allocate
stuff. This should help when combining multiple models together.
"""
struct AddModel{T,T1<:AbstractModel{T},T2<:AbstractModel{T}} <: AbstractModel{T}
    m1::T1
    m2::T2
end
Base.:+(m1::T1, m2::T2) where {T1<:AbstractModel, T2<:AbstractModel} = AddModel(m1, m2)
add(m1::M1, m2::M2) where {M1<:AbstractModel, M2<:AbstractModel} = AddModel(m1, m2)


"""
    $(SIGNATURES)
Returns the components for a composite model. This
will return a Tuple with all the models you have constructed.
"""
components(m::AbstractModel) = m
components(m::AddModel{M1,M2}) where
    {M1<:AbstractModel, M2<:AbstractModel} = (components(m.m1)..., components(m.m2)...)

flux(m::AddModel) = flux(m.m1) + flux(m.m2)

function intensity(m::AddModel{T1,T2}, x, y, args...) where {T1, T2}
    return intensity(m.m1, x, y, args...) + intensity(m.m2, x, y, args...)
end


@inline function _visibility(m::AddModel{T1,T2}, u, v, args...) where {T1, T2}
    ut,vt = transformuv(m, u, v)
    return _visibility(m.m1, ut, vt, args...) + visibility(m.m2, ut, vt, args...)
end
