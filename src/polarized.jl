struct PolarizedModel{I,Q,U,V}
    mi::I
    mq::Q
    mu::U
    mv::V
end

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
    rl = p.q + 1im*p.q
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

@inline function stokesvisibility(pimg::PolarizedModel, u, v)
    i = visibility(pimg.i, u, v)
    q = visibility(pimg.q, u, v)
    u = visibility(pimg.u, u, v)
    v = visibility(pimg.v, u, v)
    return StokesVector(i, q, u, v)
end

@inline function coherencymatrix(pimg::PolarizedModel, u, v)
    si = stokesvisibility(pimg, u, v)
    return convert(::CoherencyMatrix, si)
end
