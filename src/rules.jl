function ChainRulesCore.rrule(::Type{SA}, t::Tuple) where {SA<:StructArray}
    sa = SA(t)
    pt = ProjectTo(t)
    function _structarray_tuple_pullback(Δ)
        ps = getproperty.(Ref(Δ), propertynames(Δ))
        return NoTangent(), Tangent{typeof(t)}(ps...)
    end
    return sa, _structarray_tuple_pullback
end


function ChainRulesCore.rrule(::Type{SA}, t::NamedTuple{Na}) where {Na, SA<:StructArray}
    sa = SA(t)
    pt = ProjectTo(t)
    function _structarray_tuple_pullback(Δd)
        Δ = unthunk(Δd)
        ps = getproperty.(Ref(Δ), Na)
        nps = NamedTuple{Na}(ps)
        return NoTangent(), pt(Tangent{typeof(t)}(;nps...))
    end
    return sa, _structarray_tuple_pullback
end

# function (::ChainRulesCore.ProjectTo{T})(dx::ChainRulesCore.Tangent{<:T}) where T<:ComradeBase.AxisKeys.NamedDims.NamedDimsArray
#     println(typeof(dx.data))
#     throw("HERE")
# end

# using StructArrays
function ChainRulesCore.ProjectTo(x::StructArray)
    ProjectTo{StructArray}(;eltype = eltype(x), names=propertynames(x), dims=size(x))
end

# function extract_components(len, comp, names)
#     map(names) do n
#         c = getproperty(comp, n)
#         typeof(c) <: AbstractZero && return Fill(0.0, len)
#         return c
#     end
# end

# function (project::ProjectTo{StructArray{T}})(dx::Tangent{<:StructArray}) where {T}
#     comp = dx.components
#     tcomp = extract_components(project.dims, comp, propertynames(comp))
#     ret =  StructArray{T}(tcomp)
#     return ret
# end
function (project::ProjectTo{StructArray})(dx::Tangent{<:StructArray})
    StructArray{project.eltype}(dx.components)
end

# function (project::ProjectTo{StructArray})(dx::NamedTuple)
#     @assert project.names == keys(dx) "Key mismatch"
#     return StructArray(dx)
# end


# (project::ProjectTo{StructArray{T}})(dx::Tuple) where {T} = StructArray{T}(dx)

function (project::ProjectTo{StructArray})(dx::AbstractArray)
    @assert project.eltype === eltype(dx) "The eltype of the array is not the same there is an error in a ChainRule"
    r = StructArray(dx)
    return r
end

function (project::ProjectTo{StructArray})(dx::StructArray)
    @assert project.eltype === eltype(dx) "The eltype of the array is not the same there is an error in a ChainRule"
    return dx
end

function (project::ProjectTo{StructArray})(dx)
    @assert project.eltype === eltype(dx) "The eltype of the array is not the same there is an error in a ChainRule"
    return dx
end

function (project::ProjectTo{StructArray})(dx::ChainRulesCore.AbstractZero)
    # @assert project.eltype === eltype(dx) "The eltype of the array is not the same there is an error in a ChainRule"
    return ZeroTangent()
end




# function (project::ProjectTo{StructArray})(dx::AbstractArray) where {T}
#     # Extract the properties
#     println(typeof(dx))
#     comps = map(p->getproperty.(dx, Ref(p)), project.names)
#     StructArray{T}(comps)
# end


# function (project::ProjectTo{StructArray})(dx::AbstractArray{<:Tangent})
#     return StructArray{project.eltype}(map(p->getproperty.(dx, p), propertynames(names)))
# end


# function (project::ProjectTo{StructArray{T}})(dx::AbstractZero) where {T}
#     return dx
# end





# const AxisKeys = ComradeBase.AxisKeys
# function (project::ChainRulesCore.ProjectTo{AxisKeys.KeyedArray})(dx::ChainRulesCore.Tangent{<:AxisKeys.KeyedArray})
#     ZeroTangent()
# end

function ChainRulesCore.rrule(::Type{IntensityMap}, data::AbstractArray, keys...)
    img = IntensityMap(data, keys...)
    pd = ProjectTo(data)
    function _IntensityMap_pullback(Δ)
        return (NoTangent(), @thunk(pd(Δ)), map(i->NoTangent(), keys)...)
    end
    return img, _IntensityMap_pullback
end


function ChainRulesCore.rrule(::Type{ContinuousImage}, data::AbstractArray, pulse::AbstractModel)
    img = ContinuousImage(data, pulse)
    pd = ProjectTo(data)
    function _ContinuousImage_pullback(Δ)
        return (NoTangent(), @thunk(pd(Δ.img)), NoTangent())
    end
    return img, _ContinuousImage_pullback
end


getm(m::AbstractModel) = m
getm(m::Tuple) = m[2]

function ChainRulesCore.rrule(::typeof(visibilities_analytic), m::AbstractModel, u::AbstractArray, v::AbstractArray, t::AbstractArray, f::AbstractArray)
    vis = visibilities_analytic(m, u, v, t, f)
    function _composite_visibilities_analytic_pullback(Δ)
        du = zero(u)
        dv = zero(v)
        df = zero(f)
        dt = zero(t)

        dvis = zero(vis)
        dvis .= unthunk(Δ)
        rvis = zero(vis)
        d = autodiff(Reverse, visibilities_analytic!, Const, Duplicated(rvis, dvis), Active(m), Duplicated(u, du), Duplicated(v, dv), Duplicated(t, dt), Duplicated(f, df))
        dm = getm(d[1])
        tm = __extract_tangent(dm)
        return NoTangent(), tm, du, dv, df, dt
    end

    return vis, _composite_visibilities_analytic_pullback
end
