function ChainRulesCore.rrule(::Type{SA}, t::Tuple) where {SA<:StructArray}
    sa = SA(t)
    function _structarray_tuple_pullback(Δ)
        ps = getproperty.(Ref(Δ), propertynames(Δ))
        return NoTangent(), Tangent{typeof(t)}(map((p,tt)->p.*tt, ps, t)...)
    end
    return sa, _structarray_tuple_pullback
end


function ChainRulesCore.rrule(::Type{SA}, t::NamedTuple{Na}) where {Na, SA<:StructArray}
    sa = SA(t)
    function _structarray_tuple_pullback(Δ)
        ps = getproperty.(Ref(Δ), Na)
        ts = getproperty.(Ref(t), Na)
        nps = NamedTuple{Na}(map((p,tt)->p.*tt, ps, ts))
        return NoTangent(), Tangent{typeof(t)}(;nps...)
    end
    return sa, _structarray_tuple_pullback
end

const AxisKeys = ComradeBase.AxisKeys
function (project::ChainRulesCore.ProjectTo{AxisKeys.KeyedArray})(dx::ChainRulesCore.Tangent{<:AxisKeys.KeyedArray})
    ZeroTangent()
end

function ChainRulesCore.rrule(::Type{IntensityMap}, data::AbstractArray, keys)
    img = IntensityMap(data, keys)
    function _IntensityMap_pullback(Δ)
        return NoTangent(), parent(parent(Δ)), NoTangent()
    end
    return img, _IntensityMap_pullback
end
