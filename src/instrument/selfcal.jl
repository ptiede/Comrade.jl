export selfcal_phase

# The instrument parameter of the residual phase about the self-calibration solution.
const SELFCAL_RESIDUAL = :selfcal_residual
# Squarings of the station matrix: φ̂ is the phase of a column of A^(2^SELFCAL_SQUARINGS).
const SELFCAL_SQUARINGS = 6

struct ObservedSelfCalInstrumentModel{M <: ObservedInstrumentModel, W, G, L} <: AbstractInstrumentModel
    model::M
    phase::Symbol
    wdata::W
    weight::G
    layout::L
    squarings::Int
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
    return ObservedSelfCalInstrumentModel(intobs, SELFCAL_RESIDUAL, wdata, weight, layout, SELFCAL_SQUARINGS), prior
end

# Index maps from baselines to the station matrices `A[a, b, g]` of every group `g`, a set of
# entries of the residual at one time and frequency that baselines link, and from each entry
# to its group's reference: the entry with the largest summed weight `rowweight` over the
# baselines it is on. `pick[a, b, g]` indexes `[0; x; conj(x)]` for a per-baseline `x`: `1 + k`
# when baseline `k` runs from station `a` to `b`, `1 + K + k` when it runs from `b` to `a`, and
# `1` when no baseline joins them.
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
    refentry = [ref[ge[e]] for e in eachindex(ge)]
    column = [L[ne[e], ne[refentry[e]], ge[e]] for e in eachindex(ge)]
    eye = reshape([a == b ? 1.0 : 0.0 for a in 1:N, b in 1:N], N, N, 1)
    return (; pick, eye, column, refentry, i1 = collect(Int, i1), i2 = collect(Int, i2))
end

_nproducts(::AbstractVector{<:Number}) = 1
_nproducts(v::AbstractVector{<:StaticArraysCore.StaticMatrix}) = length(eltype(v))
_product(v::Number, ::Int) = v
_product(v::StaticArraysCore.StaticMatrix, p::Int) = v[p]

# Data and weights per correlation product; missing data get zero weight, as in the likelihood.
function _selfcal_weights(table)
    V = measurement(table)
    σ = noise(table)
    np = _nproducts(V)
    missingdata(v, s) = isnan(v) || isnan(s)
    wdata = ntuple(np) do p
        map(V, σ) do v, s
            vp, sp = _product(v, p), _product(s, p)
            missingdata(vp, sp) ? zero(complex(float(real(vp)))) : complex(vp) / sp^2
        end
    end
    weight = ntuple(np) do p
        map(V, σ) do v, s
            vp, sp = _product(v, p), _product(s, p)
            missingdata(vp, sp) ? zero(float(sp)) : inv(sp^2)
        end
    end
    return wdata, weight
end

_products(v::AbstractArray{<:Number}) = (v,)
_products(v::StructArray{<:StaticArraysCore.StaticMatrix}) = Tuple(StructArrays.components(v))
_products(v::AbstractArray{<:StaticArraysCore.StaticMatrix}) = _products(StructArray(v))
_withproducts(::AbstractArray{<:Number}, ps) = only(ps)
function _withproducts(v::AbstractArray{<:StaticArraysCore.StaticMatrix}, ps)
    T = StaticArraysCore.similar_type(eltype(v), eltype(first(ps)))
    return StructArray{T}(ps)
end

_squared(A) = dropdims(sum(reshape(A, size(A, 1), size(A, 2), 1, :) .* reshape(A, 1, size(A, 1), size(A, 2), :); dims = 2); dims = 2)
_unitnorm(A) = A ./ sqrt.(sum(abs2, A; dims = (1, 2)))

# φ̂ at every entry of the phase term, given the model visibilities without the station phase.
function _selfcal_solution(vt, J::ObservedSelfCalInstrumentModel)
    (; pick, eye, column) = J.layout
    m = _products(vt)
    c = sum(map((w, mp) -> w .* conj.(mp), J.wdata, m))
    g = sum(map((w, mp) -> w .* abs2.(mp), J.weight, m))
    cx = vcat(zero(c[1:1]), c, conj.(c))
    gx = vcat(zero(g[1:1]), g, g)
    d = sum(gx[pick]; dims = 2)
    A = cx[pick] .+ d .* eye
    for _ in 1:J.squarings
        A = _squared(_unitnorm(A))
    end
    z = A[column]
    return atan.(imag.(z), real.(z))
end

@inline function apply_instrument(vis, J::ObservedSelfCalInstrumentModel, x)
    vt = apply_instrument(vis, J.model, x)
    φ = _selfcal_solution(vt, J) .+ parent(siteparams(getproperty(x.instrument, J.phase)))
    Δ = φ[J.layout.i1] .- φ[J.layout.i2]
    r = complex.(cos.(Δ), sin.(Δ))
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
