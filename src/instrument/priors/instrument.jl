struct InstrumentPrior{A, N} <: AbstractInstrumentPrior
    site_dists::N
    function InstrumentPrior(array, site_dists::NamedTuple{N, T}) where {N, T<:NTuple{<:Any, <:InstrumentPrior}}
        dists = map(x->build_site_dist(array, x), site_dists)
        return new{typeof(array), typeof(dists)}(dist)
    end
end
