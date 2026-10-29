export selfcal_phase

# The instrument parameter of the residual phase about the self-calibration solution.
const SELFCAL_RESIDUAL = :selfcal_residual
# Squarings of the station matrix: φ̂ is the phase of a column of A^(2^SELFCAL_SQUARINGS).
const SELFCAL_SQUARINGS = 6

struct ObservedSelfCalInstrumentModel{M <: ObservedInstrumentModel, W, G, L, K, E} <: AbstractInstrumentModel
    model::M
    phase::Symbol
    wdata::W
    weight::G
    layout::L
    squarings::Int
    # the measured values, by record then channel
    keep::K
    # the residual's entries at each value's first and second station
    entry1::E
    entry2::E
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
    return InstrumentModel(model.jones, prior, model.refbasis, model.gaugefix, false, model.feedrotation_corrected)
end

# A self-calibrated model is solved against the one table it is fit to.
_checkdata(::AbstractInstrumentModel, data) = nothing
function _checkdata(m::InstrumentModel, data)
    m.selfcal || return nothing
    vis = filter(d -> d isa VisTable && !(eltype(measurement(d)) <: Real), collect(data))
    length(data) == 1 && length(vis) == 1 || throw(
        ArgumentError(
            "a self-calibrated instrument model needs exactly one table of visibilities or coherencies, " *
                "got $(length(data)) table(s), $(length(vis)) of them complex"
        )
    )
    return nothing
end

function _set_selfcal(m::InstrumentModel, t::VisTable)
    keep = _measuredvalues(t)
    wdata, weight = _selfcal_weights(t, keep)
    # The entries of the residual and the baselines between them do not depend on which are
    # fixed, so a first build without references gives the groups the references come from.
    intfree, priorfree = set_array(_with_residual(m, NoReference()), t)
    smap = sitelookup(getproperty(priorfree, SELFCAL_RESIDUAL))
    e1, e2 = _entrymaps(intfree, SELFCAL_RESIDUAL)
    rowweight = sum(map(w -> collect(w), weight))
    layout = _selfcal_layout(smap, e1[keep], e2[keep], rowweight)
    refs = unique(layout.refentry)
    intobs, prior = set_array(_with_residual(m, EntryReference(refs, 0.0)), t)
    sort(fixed_indices(getproperty(prior, SELFCAL_RESIDUAL))) == sort(refs) ||
        error("the residual phase prior did not fix exactly the reference entries")
    entry1, entry2 = _entrymaps(intobs, SELFCAL_RESIDUAL)
    return ObservedSelfCalInstrumentModel(
        intobs, SELFCAL_RESIDUAL, wdata, weight, layout, SELFCAL_SQUARINGS, keep, entry1, entry2
    ), prior
end

# The entries of the term `name` at each value's first and second station.
function _entrymaps(J::ObservedInstrumentModel, name)
    L = layout(J)
    e = getproperty(L.entries, name)
    return e[L.cell1], e[L.cell2]
end

# Index maps from baselines to the station matrices of every group, a set of entries of the
# residual at one time and frequency that baselines link, and from each entry to its group's
# reference: the entry with the largest summed weight `rowweight` over the baselines it is on.
# The groups with `n` stations form one stack, whose `pick[a, b, j]` indexes `[0; x; conj(x)]`
# for a per-baseline `x`: `1 + k` when baseline `k` runs from station `a` to `b` of the stack's
# `j`-th group, `1 + K + k` when it runs from `b` to `a`, and `1` when no baseline joins them.
# An entry's phase is that of `A[column]` for its stack's matrices `A`, and `order` sorts the
# entries of the stacks, concatenated, into the order of the residual.
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
    G = length(group_of)
    groupsites = [sort(unique(st[e] for e in eachindex(ge) if ge[e] == g)) for g in 1:G]
    ne = [findfirst(==(st[e]), groupsites[ge[e]]) for e in eachindex(st)]
    nstations = map(length, groupsites)

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
    refentry = [ref[ge[e]] for e in eachindex(ge)]

    sizes = sort(unique(nstations))
    members = [findall(==(n), nstations) for n in sizes]
    stack = [searchsortedfirst(sizes, n) for n in nstations]
    slot = zeros(Int, G)
    for m in members
        slot[m] .= eachindex(m)
    end

    K = length(i1)
    picks = [ones(Int, n, n, length(m)) for (n, m) in zip(sizes, members)]
    for k in eachindex(i1, i2)
        e1, e2 = i1[k], i2[k]
        g, a, b = ge[e1], ne[e1], ne[e2]
        pick = picks[stack[g]]
        pick[a, b, slot[g]] == 1 ||
            throw(ArgumentError("baseline $k repeats the station pair of $(describe(e1)) and $(describe(e2))"))
        pick[a, b, slot[g]] = 1 + k
        pick[b, a, slot[g]] = 1 + K + k
    end

    entries = [findall(e -> stack[ge[e]] == j, eachindex(ge)) for j in eachindex(sizes)]
    stacks = map(eachindex(sizes)) do j
        L = LinearIndices(picks[j])
        column = [L[ne[e], ne[refentry[e]], slot[ge[e]]] for e in entries[j]]
        n = sizes[j]
        eye = reshape([a == b ? 1.0 : 0.0 for a in 1:n, b in 1:n], n, n, 1)
        (; pick = picks[j], eye, column)
    end
    order = invperm(reduce(vcat, entries))
    return (; stacks = Tuple(stacks), order, refentry)
end

# Data and weights per correlation product at the measured values `keep`.
function _selfcal_weights(t::VisTable, keep)
    V, w = _dense(measurement(t)), _dense(weight(t))
    Vs, ws = _feedproducts(V), _feedproducts(w)
    wdata = map((v, wk) -> wk[keep] .* v[keep], Vs, ws)
    return wdata, map(wk -> wk[keep], ws)
end

# The correlation products of dense values: itself, or the four feed slabs `(…, 2, 2)`.
_feedproducts(x::AbstractArray{<:Number, 2}) = (x,)
_feedproducts(x::AbstractArray{<:Number, 4}) = ntuple(k -> x[:, :, mod1(k, 2), cld(k, 2)], Val(4))

_products(v::AbstractArray{<:Number}) = (v,)
_products(v::FieldDimArray{<:StaticArraysCore.StaticMatrix}) = ntuple(k -> fieldview(v, k), Val(4))

# φ̂ at every entry of the phase term, given the model visibilities without the station phase.
function _selfcal_solution(vt, J::ObservedSelfCalInstrumentModel)
    m = map(p -> p[J.keep], _products(vt))
    c = sum(map((w, mp) -> w .* conj.(mp), J.wdata, m))
    g = sum(map((w, mp) -> w .* abs2.(mp), J.weight, m))
    cx = vcat(zero(c[1:1]), c, conj.(c))
    gx = vcat(zero(g[1:1]), g, g)
    z = map(J.layout.stacks) do st
        A = cx[st.pick] .+ sum(gx[st.pick]; dims = 2) .* st.eye
        return _normsquares(A, J.squarings)[st.column]
    end
    zo = reduce(vcat, z)[J.layout.order]
    return atan.(imag.(zo), real.(zo))
end

# `squarings` times `A ← (A / ‖A‖)²` on each matrix of the stack `A` of Hermitian matrices,
# with `‖⋅‖` the Frobenius norm.
function _normsquares(A, squarings)
    for _ in 1:squarings
        A = _squared(_unitnorm(A))
    end
    return A
end

_unitnorm(A) = A ./ sqrt.(sum(abs2, A; dims = (1, 2)))
_squared(A) = _matmul(A, A)
_matmul(A, B) = dropdims(sum(reshape(A, size(A, 1), size(A, 2), 1, :) .* reshape(B, 1, size(B, 1), size(B, 2), :); dims = 2); dims = 2)
_adjoint(A) = conj.(permutedims(A, (2, 1, 3)))

# Adjoints are `x̄ = ∂L/∂re(x) + i ∂L/∂im(x)`. The adjoint of `_normsquares` at the powers
# `[A_0, …, A_S]` for the output adjoint `Ȳ` holds for Hermitian perturbations of `A_0`, so
# each step's adjoint is projected onto the Hermitian matrices. For `B = s A²`,
# `s = 1 / ‖A‖²`: `Ā = s (B̄ A + A B̄) - 2 s re⟨B̄, B⟩ A`.
function _normsquares_adjoint(Ȳ, powers)
    B̄ = (Ȳ .+ _adjoint(Ȳ)) ./ 2
    for i in lastindex(powers):-1:(firstindex(powers) + 1)
        A, B = powers[i - 1], powers[i]
        s = inv.(sum(abs2, A; dims = (1, 2)))
        r = real.(sum(conj.(B̄) .* B; dims = (1, 2)))
        X = _matmul(B̄, A)
        B̄ = s .* (X .+ _adjoint(X)) .- 2 .* s .* r .* A
    end
    return B̄
end

function EnzymeRules.augmented_primal(
        config::EnzymeRules.RevConfigWidth, ::Const{typeof(_normsquares)}, ::Type{<:Annotation},
        A::Annotation, squarings::Const
    )
    powers = [EnzymeRules.overwritten(config)[2] ? copy(A.val) : A.val]
    for _ in 1:squarings.val
        push!(powers, _squared(_unitnorm(last(powers))))
    end
    Y = last(powers)
    W = EnzymeRules.width(config)
    shadow = EnzymeRules.needs_shadow(config) ? (W == 1 ? zero(Y) : ntuple(_ -> zero(Y), Val(W))) : nothing
    primal = EnzymeRules.needs_primal(config) ? Y : nothing
    return EnzymeRules.AugmentedReturn(primal, shadow, (powers, shadow))
end

function EnzymeRules.reverse(
        config::EnzymeRules.RevConfigWidth, ::Const{typeof(_normsquares)}, ::Type{<:Annotation}, tape,
        A::Annotation, squarings::Const
    )
    powers, shadow = tape
    W = EnzymeRules.width(config)
    for w in 1:W
        Ā = _shadow(A, w, W)
        (isnothing(Ā) || isnothing(shadow)) && continue
        Ā .+= _normsquares_adjoint(W == 1 ? shadow : shadow[w], powers)
    end
    return (nothing, nothing)
end

# The `w`-th of `W` shadows of `x`, or `nothing` when `x` is inactive.
_shadow(::Const, w, W) = nothing
function _shadow(x::Annotation, w, W)
    dx = W == 1 ? x.dval : x.dval[w]
    return dx === x.val ? nothing : dx
end

@inline function apply_instrument(vis, J::ObservedSelfCalInstrumentModel, x)
    vt = apply_instrument(vis, J.model, x)
    φ = _selfcal_solution(vt, J) .+ parent(siteparams(getproperty(x.instrument, J.phase)))
    Δ = φ[J.entry1] .- φ[J.entry2]
    return _scaled(vt, complex.(cos.(Δ), sin.(Δ)))
end

_scaled(v::AbstractArray{<:Number}, r) = v .* r
_scaled(v::FieldDimArray{<:StaticArraysCore.StaticMatrix}, r) = FieldDimArray{SMatrix{2, 2}}(parent(v) .* r)

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
