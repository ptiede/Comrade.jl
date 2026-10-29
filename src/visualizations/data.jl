# Unmeasured values and hands are `NaN`.
_rowvalue(m, σ) = isfinite(σ) ? (m, σ) : (oftype(m, NaN), oftype(σ, NaN))
function _rowvalue(m::StaticArraysCore.StaticMatrix, σ::StaticArraysCore.StaticMatrix)
    v = map(_rowvalue, m, σ)
    return map(first, v), map(last, v)
end

export datarows, residual

"""
    datarows(t)

The measured values of the [`VisTable`](@ref) or [`ClosureTable`](@ref) `t` as a table with
one row each. For visibilities the columns are `U`, `V` (λ), `uvdist`, `Ti`, `Fr`, `sites`,
`measurement` and `noise`; the unmeasured correlations of a `Correlations()` table are `NaN`. For closures they are
`uvdist` (the longest baseline of the closure, λ), `Ti`, `Fr`, `stations`, `measurement` and
`noise`.
"""
function datarows(t::VisTable)
    keep = _measuredvalues(t)
    d = domain(t)
    U, V = d.U[keep], d.V[keep]
    (; Ti, Fr, sites) = _measuredpoints(t)
    mv = map(_rowvalue, baseimage(measurement(t))[keep], baseimage(noise(t))[keep])
    return StructArray((;
        U, V, uvdist = hypot.(U, V), Ti, Fr, sites, measurement = map(first, mv), noise = map(last, mv),
    ))
end

function datarows(ct::ClosureTable)
    vis = getfield(ct, :vis)
    keep = _measuredvalues(vis)
    uvd = hypot.(domain(vis).U[keep], domain(vis).V[keep])
    Dt = sparse(getfield(ct, :design)')
    longest = [maximum(uvd[rowvals(Dt)[nzrange(Dt, r)]]) for r in axes(Dt, 2)]
    d = domain(ct)
    return StructArray((;
        uvdist = longest, Ti = d.Ti, Fr = d.Fr, stations = d.stations,
        measurement = Vector(baseimage(measurement(ct))), noise = Vector(baseimage(noise(ct))),
    ))
end

"""
    residual(post::AbstractVLBIPosterior, p)

A figure of the normalized residuals of every data product of the posterior `post` at the
parameters `p`, one row per product, against the baseline length. Needs a Makie backend.
"""
function residual end

