export selfcal_phase

# The instrument parameter of the residual phase about the self-calibration solution.
const SELFCAL_RESIDUAL = :selfcal_residual
# Power iterations on the station matrix: φ̂ is the phase of A^SELFCAL_ITERATIONS e_ref.
const SELFCAL_ITERATIONS = 32

struct ObservedSelfCalInstrumentModel{M <: ObservedInstrumentModel, W, G, L} <: AbstractInstrumentModel
    model::M
    phase::Symbol
    wdata::W
    weight::G
    layout::L
    iterations::Int
end

EnzymeRules.inactive_type(::Type{<:ObservedSelfCalInstrumentModel}) = true

function Base.show(io::IO, mime::MIME"text/plain", m::ObservedSelfCalInstrumentModel)
    printstyled(io, "ObservedSelfCalInstrumentModel"; bold = true, color = :light_cyan)
    println(io, " (station phase = self-calibration + $(m.phase))")
    return show(io, mime, m.model)
end

# `model` with the residual phase added, referenced by `refant`.
function _with_residual(model::InstrumentModel, refant)
    # Uniform on the circle to within 3e-9, with one latent coordinate per entry, no singular
    # point, and no seam at ±π.
    residual = ArrayPrior(IIDSitePrior(IntegSeg(), WrappedNormal(0.0, 2π)); refant)
    prior = merge(model.prior, NamedTuple{(SELFCAL_RESIDUAL,)}((residual,)))
    return InstrumentModel(model.jones, prior, model.refbasis, model.gaugefix, false)
end

function set_array(m::InstrumentModel, array::AbstractArrayConfiguration, data)
    m.selfcal || return set_array(m, array)
    tables = filter(d -> d isa EHTObservationTable{<:Union{EHTCoherencyDatum, EHTVisibilityDatum}}, collect(data))
    length(tables) == 1 || throw(
        ArgumentError(
            "a self-calibrated instrument model needs exactly one table of visibilities or coherencies, " *
                "got $(length(tables))"
        )
    )
    table = only(tables)
    length(table) == length(array) || throw(
        DimensionMismatch("the data have $(length(table)) rows but the array has $(length(array))")
    )
    wdata, weight = _selfcal_weights(table)
    # The entries of the residual and the baselines between them do not depend on which are
    # fixed, so a first build without references gives the groups the references come from.
    intfree, priorfree = set_array(_with_residual(m, NoReference()), array)
    smap = sitelookup(getproperty(priorfree, SELFCAL_RESIDUAL))
    bl = getproperty(intfree.bsitelookup, SELFCAL_RESIDUAL)
    rowweight = sum(map(w -> collect(w), weight))
    layout = _selfcal_layout(smap, bl.indices_1, bl.indices_2, rowweight)
    refs = unique(layout.refentry)
    intobs, prior = set_array(_with_residual(m, EntryReference(refs, 0.0)), array)
    sort(fixed_indices(getproperty(prior, SELFCAL_RESIDUAL))) == sort(refs) ||
        error("the residual phase prior did not fix exactly the reference entries")
    return ObservedSelfCalInstrumentModel(intobs, SELFCAL_RESIDUAL, wdata, weight, layout, SELFCAL_ITERATIONS), prior
end

# Index maps from baselines to the station matrices `A[a, b, g]` of every group `g`, a set of
# entries of the residual at one time and frequency that baselines link, and from each entry
# to its group's reference: the entry with the largest summed weight `rowweight` over the
# baselines it is on. `pick[a, b, g]` indexes `[0; x; conj(x)]` for a per-baseline `x`: `1 + k`
# when baseline `k` runs from station `a` to `b`, `1 + K + k` when it runs from `b` to `a`, and
# `1` when no baseline joins them. `start[a, g]` indexes `A[a, r, g]` for the reference station
# `r` of group `g`, `entry[e]` indexes the `N × G` station vectors at entry `e`, and
# `refslot[e]` at the reference entry of `e`'s group.
function _selfcal_layout(smap::SiteLookup, i1::AbstractVector{<:Integer}, i2::AbstractVector{<:Integer}, rowweight)
    st = sites(smap)
    stamps = collect(zip(times(smap), frequencies(smap)))
    describe(e) = "site = $(st[e]), t = $(_center(times(smap)[e]))"

    root = collect(eachindex(st))
    function findroot(e)
        while root[e] != e
            root[e] = root[root[e]]
            e = root[e]
        end
        return e
    end
    for k in eachindex(i1, i2)
        e1, e2 = i1[k], i2[k]
        stamps[e1] == stamps[e2] || throw(
            ArgumentError("baseline $k joins $(describe(e1)) and $(describe(e2)) at different times")
        )
        root[findroot(e1)] = findroot(e2)
    end
    roots = map(findroot, eachindex(st))
    group_of = Dict(r => g for (g, r) in pairs(unique(roots)))
    ge = [group_of[r] for r in roots]
    stations = sort(unique(st))
    ne = [findfirst(==(s), stations) for s in st]
    N, G = length(stations), length(group_of)

    entryweight = zeros(length(st))
    for k in eachindex(i1, i2, rowweight)
        entryweight[i1[k]] += rowweight[k]
        entryweight[i2[k]] += rowweight[k]
    end
    ref = zeros(Int, G)
    for e in eachindex(ge)
        g = ge[e]
        (ref[g] == 0 || entryweight[e] > entryweight[ref[g]]) && (ref[g] = e)
    end

    K = length(i1)
    pick = ones(Int, N, N, G)
    for k in eachindex(i1, i2)
        e1, e2 = i1[k], i2[k]
        g, a, b = ge[e1], ne[e1], ne[e2]
        pick[a, b, g] == 1 ||
            throw(ArgumentError("baseline $k repeats the station pair of $(describe(e1)) and $(describe(e2))"))
        pick[a, b, g] = 1 + k
        pick[b, a, g] = 1 + K + k
    end

    L = LinearIndices((N, N, G))
    start = [L[a, ne[ref[g]], g] for a in 1:N, g in 1:G]
    entry = [LinearIndices((N, G))[ne[e], ge[e]] for e in eachindex(ge)]
    refentry = [ref[ge[e]] for e in eachindex(ge)]
    eye = reshape([a == b ? 1.0 : 0.0 for a in 1:N, b in 1:N], N, N, 1)
    return (; pick, eye, start, entry, refslot = entry[refentry], refentry, i1 = collect(Int, i1), i2 = collect(Int, i2))
end

_products(v::AbstractArray{<:Number}) = (v,)
_products(v::StructArray{<:StaticArraysCore.StaticMatrix}) = Tuple(StructArrays.components(v))
_products(v::AbstractArray{<:StaticArraysCore.StaticMatrix}) = _products(StructArray(v))
_withproducts(::AbstractArray{<:Number}, ps) = only(ps)
function _withproducts(v::AbstractArray{<:StaticArraysCore.StaticMatrix}, ps)
    T = StaticArraysCore.similar_type(eltype(v), eltype(first(ps)))
    return StructArray{T}(ps)
end

# Data and weights per correlation product; missing data get zero weight, as in the likelihood.
function _selfcal_weights(table)
    V, σ = _products(measurement(table)), _products(noise(table))
    observed(v, s) = !(isnan(v) || isnan(s))
    wdata = map(V, σ) do Vp, σp
        map((v, s) -> observed(v, s) ? complex(v) / s^2 : zero(complex(float(real(v)))), Vp, σp)
    end
    weight = map(V, σ) do Vp, σp
        map((v, s) -> observed(v, s) ? inv(s^2) : zero(float(s)), Vp, σp)
    end
    return wdata, weight
end

_matvec(A, x) = dropdims(sum(A .* reshape(x, 1, size(x)...); dims = 2); dims = 2)
_unitnorm(x) = x ./ sqrt.(sum(abs2, x; dims = 1))

# φ̂ at every entry of the phase term, given the model visibilities without the station phase.
function _selfcal_solution(vt, J::ObservedSelfCalInstrumentModel)
    (; pick, eye, start, entry, refslot) = J.layout
    m = _products(vt)
    c = sum(map((w, mp) -> w .* conj.(mp), J.wdata, m))
    g = sum(map((w, mp) -> w .* abs2.(mp), J.weight, m))
    cx = vcat(zero(c[1:1]), c, conj.(c))
    gx = vcat(zero(g[1:1]), g, g)
    d = sum(gx[pick]; dims = 2)
    A = cx[pick] .+ d .* eye
    x = A[start]
    for _ in 2:J.iterations
        x = _matvec(A, _unitnorm(x))
    end
    z = x[entry] .* conj.(x[refslot])
    return atan.(imag.(z), real.(z))
end

@inline function apply_instrument(vis, J::ObservedSelfCalInstrumentModel, x)
    vt = apply_instrument(vis, J.model, x)
    φ = _selfcal_solution(vt, J) .+ parent(siteparams(getproperty(x.instrument, J.phase)))
    Δ = φ[J.layout.i1] .- φ[J.layout.i2]
    r = cis.(Δ)
    return _withproducts(vt, map(mp -> mp .* r, _products(vt)))
end

"""
    selfcal_phase(post::VLBIPosterior, θ)

The station phase `φ = φ̂ + residual` at `θ` of a posterior whose instrument model
self-calibrates the station phase (see [`InstrumentModel`](@ref)), as a `SiteArray`.
"""
function selfcal_phase(post, θ)
    J = instrumentmodel(post)
    J isa ObservedSelfCalInstrumentModel ||
        throw(ArgumentError("the posterior's instrument model does not self-calibrate the station phase"))
    _, vis = idealmaps(datatype(typeof(post)), skymodel(post), θ)
    vt = apply_instrument(vis, J.model, θ)
    ε = siteparams(getproperty(θ.instrument, J.phase))
    return SiteArray(_selfcal_solution(vt, J) .+ parent(ε), times(ε), frequencies(ε), sites(ε))
end
