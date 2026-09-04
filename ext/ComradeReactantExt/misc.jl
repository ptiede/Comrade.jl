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
#
# `y[inds] .= vals` lowers to a `stablehlo.scatter` — the fast, vectorized ("raised")
# path — for length > 1, but its single-element case hits a Reactant bug that tries to
# assign a 1-element traced array into a scalar slot (`Float64(::TracedRArray)`), e.g. a
# TrackSeg offset referenced to one antenna has exactly one fixed index. `inds` is a
# concrete host vector, so the length branches are resolved at trace time (no traced
# branch): empty is a no-op and only the length-1 case falls back to a scalar store.
# See EnzymeAD/Reactant.jl#2960.
@inline function Comrade.scatter_values!(y::Reactant.AnyTracedRArray, inds, vals)
    isempty(inds) && return y
    Reactant.@allowscalar begin
        if length(inds) == 1
            y[inds[begin]] = vals[begin]
        else
            y[inds] .= vals
        end
    end
    return y
end

# Log-depth (Hillis-Steele) scan for the recurrence in the generic method. Affine
# maps compose associatively — `(a₂,c₂) ∘ (a₁,c₁) = (a₂a₁, a₂c₁ + c₂)` — so after
# the round with shift `s` entries `1:s` are final, and one slice-multiply-add folds
# the partial scan `s` places into the tail, doubling the resolved prefix. The
# `while` runs on the HOST at trace time and unrolls into ⌈log₂ n⌉ rounds of dense
# ops (~13 rounds for IntegSeg-length chains, hard-capped at 64) — it exists because
# Reactant's loop raising never vectorizes a loop-carried recurrence, so the
# sequential form would serialize into a `stablehlo.while` in every logdensity and
# gradient evaluation. Composition is plain mul/add, so the exact zeros the coloring
# places in `a` at chain-opening and fixed-left points reset the carry across the
# batched chain concatenation with exact gradients. Mutates its arguments
# (trace-local staging temps at both call sites).
#
# NB: the prettier work-efficient formulations are blocked upstream — the recursive
# odd-even (Blelloch) scan needs differently-sized subslices of one array, whose
# cotangents Enzyme-MLIR currently mis-accumulates (`stablehlo.add` shape mismatch),
# and `cumprod`-based closed forms are unstable with a broken Enzyme gradient. If
# Reactant ships an associative-scan primitive, delete this method — the generic
# sequential definition then serves both backends.
function Comrade.affine_scan!(a::Reactant.AnyTracedRArray, c::Reactant.AnyTracedRArray)
    n = length(c)
    s = 1
    while s < n
        at = a[(s + 1):n]
        c[(s + 1):n] = at .* c[1:(n - s)] .+ c[(s + 1):n]
        a[(s + 1):n] = at .* a[1:(n - s)]
        s <<= 1
    end
    return c
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

# Batched chain groups are worth their per-point parameter gather only where a per-chain
# walker would cost a kernel launch per chain — i.e. exactly when the parameter vector is
# traced. On the host the per-chain form wins (see `Comrade._walk_units`).
@inline function Comrade._walk_units(d::Comrade.GaussMarkovChainDist, ::Reactant.AnyTracedRArray)
    return Comrade.chaingroups(d)
end
