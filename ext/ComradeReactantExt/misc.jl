@inline function Comrade._apply_instrument!(
        vout::Union{Reactant.AnyTracedRArray, AbstractArray{<:SMatrix{2, 2, <:Reactant.TracedRNumber}}},
        vis,
        J::Comrade.ObservedInstrumentModel,
        xint
    )
    Reactant.@allowscalar @trace track_numbers = false for i in eachindex(vis)
        vout[i] = Comrade.apply_jones(vis[i], i, J, xint)
    end
    return vout
end

# Inside the traced loop the datum index is a TracedRNumber: the Ti/Fr reads are gathers
# from constant host vectors (which trace), but a Symbol cannot be gathered, so the
# JonesPoint's site is `nothing` under compilation (documented on `JonesPoint`).
@inline function Comrade._sitepoint(::Val{N}, meta, i::Reactant.TracedRNumber) where {N}
    return Comrade.JonesPoint{N}(
        Comrade.rgetindex(meta.Ti, i), Comrade.rgetindex(meta.Fr, i), nothing, i
    )
end

function StructArrays.createinstance(::Type{<:StokesParams}, args...)
    return StokesParams(args...)
end

# The indexed assignment lowers to `stablehlo.scatter`, but its bounds check iterates the
# index vector with scalar `getindex`, which is disallowed while tracing. Allowing scalar
# indexing here lets the trace through. See https://github.com/EnzymeAD/Reactant.jl/issues/2960.
@inline function Comrade.fill_partially_fixed!(
        yfv::Reactant.AnyTracedRArray, variate_index, fixed_index, y, fixed_values
    )
    Reactant.@allowscalar begin
        yfv[variate_index] = y
        yfv[fixed_index] .= fixed_values
    end
    return yfv
end

# Chained-prior gather/scatter between the flat storage and per-chain vectors: views
# with vector indices do not trace, so materialize a `stablehlo.gather`/`scatter`. The
# `@allowscalar` is for the bounds check, as in `fill_partially_fixed!` above.
@inline function Comrade.rgather(x::Reactant.AnyTracedRArray, inds)
    return Reactant.@allowscalar x[inds]
end
@inline function Comrade.rscatter!(y::Reactant.AnyTracedRArray, inds, vals)
    Reactant.@allowscalar y[inds] = vals
    return y
end

