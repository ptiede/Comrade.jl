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
        # `yfv[fixed_index] .= fixed_values` lowers to a `stablehlo.scatter` — the fast, vectorized
        # ("raised") path — for length > 1, but its single-element case hits a Reactant bug that
        # tries to assign a 1-element traced array into a scalar slot (`Float64(::TracedRArray)`),
        # e.g. a TrackSeg offset referenced to one antenna has exactly one fixed index. `fixed_index`
        # is a concrete host vector, so this length branch is resolved at trace time (no traced
        # branch) and only the length-1 case falls back to a scalar store. See EnzymeAD/Reactant.jl#2960.
        if length(fixed_index) == 1
            yfv[fixed_index[begin]] = fixed_values[begin]
        else
            yfv[fixed_index] .= fixed_values
        end
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
