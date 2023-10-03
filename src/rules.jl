
#from Lux to speed up tuple merging
function ChainRulesCore.rrule(::typeof(merge), nt1::NamedTuple{F1}, nt2::NamedTuple{F2}) where {F1, F2}
    y = merge(nt1, nt2)
    function ∇merge(dy)
        dnt1 = NamedTuple((f1 => (f1 in F2 ? NoTangent() : getproperty(dy, f1)) for f1 in F1))
        dnt2 = NamedTuple((f2 => getproperty(dy, f2) for f2 in F2))
        return (NoTangent(), dnt1, dnt2)
    end
    function ∇merge(dy::Union{NoTangent, ZeroTangent})
        return (NoTangent(), NoTangent(), NoTangent())
    end
    return y, ∇merge
end

function ChainRulesCore.rrule(::typeof(vec), x::AbstractMatrix)
    y = vec(x)
    ∇vec(dy) = (NoTangent(), reshape(dy, size(x)))
    return y, ∇vec
end

function ChainRulesCore.rrule(::typeof(collect), v::Vector)
    y = collect(v)
    ∇collect(dy) = (NoTangent(), dy)
    return y, ∇collect
end

function ChainRulesCore.rrule(::typeof(copy), x)
    ∇copy(dy) = (NoTangent(), dy)
    return copy(x), ∇copy
end
