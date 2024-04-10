"""
    $(TYPEDEF)

Holds the pairs of Jones matrices for the first and second sites of a baseline.

# Fields
$(FIELDS)
"""
struct JonesPairs{T, M1<:AbstractVector{T}, M2<:AbstractVector{T}}
    """
    Vector of jones matrices for sites 1
    """
    m1::M1
    """
    Vector of jones matrices for sites 2
    """
    m2::M2
end


struct JonesPairsStyle <: Broadcast.BroadcastStyle end
Base.BroadcastStyle(::Type{<:JonesPairs}) = JonesPairsStyle()
Base.broadcastable(x::JonesPairs) = x

function Base.similar(bc::Broadcast.Broadcasted{JonesPairsStyle}, ::Type{ElType}) where {ElType}
    m1 = similar(Vector{ElType}, length(bc))
    m2 = similar(Vector{ElType}, length(bc))
    JonesPairs(m1, m2)
end


function Base.copyto!(dest::JonesPairs, bc::Broadcast.Broadcasted{JonesPairsStyle})
    # TODO check if this is reasonable
    fbc = Broadcast.flatten(bc)
    m1 = map(Base.Fix2(getproperty, :m1), fbc.args)
    m2 = map(Base.Fix2(getproperty, :m2), fbc.args)
    _jonesmap!(fbc.f, dest.m1, dest.m2, m1, m2)
    return dest
end

@inline Base.eltype(::JonesPairs{T}) where {T} = T
@inline Base.length(j::JonesPairs) = length(j.m1)
@inline Base.size(j::JonesPairs) = (length(j),)
Base.getindex(j::JonesPairs, i::Int) = (j.m1[i], j.m2[i])
function Base.setindex!(j::JonesPairs, X, i::Int)
    j.m1[i] = X[1]
    j.m2[i] = X[2]
    return j
end
Base.IndexStyle(::Type{JonesPairs}) = IndexLinear()
Base.similar(j::JonesPairs, ::Type{S}, dims::Dims{1}) where {S} = JonesPairs(similar(j.m1, S, dims), similar(j.m2, S, dims))
Base.firstindex(j::JonesPairs) = firstindex(j.m1)
Base.lastindex(j::JonesPairs) = lastindex(j.m1)
