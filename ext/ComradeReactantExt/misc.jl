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

@inline function Comrade.site_sum(y::Reactant.AnyTracedRArray, site_map::Comrade.SiteLookup)
    yout = similar(y)
    vals = values(lookup(site_map))
    # `vals` is a Tuple of per-site index vectors (one entry per site), known at trace
    # time, so we unroll it with a plain `for` — a `@trace` loop expects a numeric range
    # and would call `step(::Tuple)`. Each `y[site]`/`yout[site] = …` is a gather/scatter
    # over the site's index vector, which traces.
    for site in vals
        ys = y[site]
        yout[site] = cumsum(ys)
    end
    return yout
end

function Comrade.branchcut(x::Reactant.TracedRNumber)
    xmod = mod(x, oftype(x, 2π))
    return ifelse(xmod > π, xmod - oftype(x, 2π), xmod)
end
