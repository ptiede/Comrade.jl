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

# Reactant's loop raising vectorizes a traced loop only when it carries a single value, so a
# matrix-valued output fills each of its component arrays with its own loop.
@inline function Comrade._apply_instrument!(
        vout::StructArray{<:SMatrix{2, 2, <:Reactant.TracedRNumber}},
        vis,
        J::Comrade.ObservedInstrumentModel,
        xint
    )
    foreach(enumerate(StructArrays.components(vout))) do (c, comp)
        Reactant.@allowscalar @trace track_numbers = false for i in eachindex(vis)
            comp[i] = Comrade.apply_jones(vis[i], i, J, xint)[c]
        end
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
# A single index fails in the broadcast form (`Float64(::TracedRArray)`), e.g. a TrackSeg
# offset referenced to one antenna, so it uses a scalar store. `inds` is a host vector, so
# the length checks resolve at trace time.
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

# Log-depth (Hillis–Steele) scan of the affine recurrence of the generic method, using
# `(a₂, c₂) ∘ (a₁, c₁) = (a₂a₁, a₂c₁ + c₂)`. The `while` loop runs at trace time and unrolls
# into ⌈log₂ n⌉ rounds of dense operations; a sequential recurrence would trace to a
# `stablehlo.while` that Reactant does not vectorize. Zeros in `a` reset the carry between
# concatenated chains, with exact gradients. Mutates both arguments.
#
# A Blelloch scan is not used: Enzyme-MLIR mis-accumulates the cotangents of its
# differently-sized subslices. `cumprod` closed forms are unstable. If Reactant provides an
# associative-scan primitive, this method can be deleted.
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

# Traced parameters walk batched chain groups: a per-chain walk would launch one kernel per
# chain. The host walks chains one at a time.
@inline function Comrade._walk_units(d::Comrade.GaussMarkovChainDist, ::Reactant.AnyTracedRArray)
    return Comrade.chaingroups(d)
end
