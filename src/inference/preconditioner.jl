# Low-rank affine preconditioning of the flat latent space the samplers work in.
#
# A diagonal mass matrix is a per-coordinate rescale and cannot align with posterior
# correlations that couple parameter blocks (station gain phases trading against image
# structure being the canonical VLBI example). The preconditioner supplies the missing
# rotation: NUTS samples z, and the model sees x = b + A z with
#
#     A = Diagonal(d) * (I + V * Diagonal(s .- 1) * V'),    V'V = I  (n × m, m ≪ n),
#
# so the m leading correlation directions of a pilot posterior are rescaled to unit
# width while the other n − m directions pass through untouched and stay the diagonal
# metric's job. A is square: the latent dimension is unchanged, and log|det A| =
# Σ log d + Σ log s is a constant.

export LowRankPreconditioner, fit_preconditioner

# Parametric so the arrays can be host `Array`s (baked as constants under Reactant)
# or `ConcreteRArray`s (traced as runtime inputs, updatable in place between compiled
# calls — the basis of recompile-free warmup refits). Validation runs only for host
# arrays; device buffers are produced from an already-validated host fit.
"""
    LowRankPreconditioner(b, d, V, s)

Affine reparameterization `x = b .+ d .* (z + V*((s .- 1) .* (V'z)))` of a sampler's
flat latent space: `b` and `d` are the per-coordinate center and scale (length `n`),
`V` an `n × m` matrix with orthonormal columns, and `s` the `m` per-direction scale
factors applied on `span(V)`.

Pass it as the `transport_method` of [`sample`](@ref); `transport_to` composes it in
front of the posterior's [`asflat`](@ref asflat) transform, so the result is an ordinary
[`TransformedVLBIPosterior`](@ref). Build it from a pilot run with
[`fit_preconditioner`](@ref).

Fits produce plain `Float64` arrays, which serialize stably into a run's
`transport.jls` and enter a Reactant-compiled program as constants. With Reactant loaded
the buffers can instead be `ConcreteRArray`s, which trace as runtime inputs and can be
updated in place between compiled calls (recompile-free warmup refits).
"""
struct LowRankPreconditioner{TB <: AbstractVector, TV_ <: AbstractMatrix, TS <: AbstractVector}
    b::TB
    d::TB
    V::TV_
    s::TS
    function LowRankPreconditioner(b::TB, d::TB, V::TV_, s::TS) where {TB, TV_, TS}
        n, m = size(V)
        (length(b) == n && length(d) == n && length(s) == m) || throw(
            DimensionMismatch(
                "b ($(length(b))), d ($(length(d))), V ($(size(V))), s ($(length(s))) " *
                    "are inconsistent: need length(b) == length(d) == size(V, 1) and " *
                    "length(s) == size(V, 2)"
            )
        )
        if b isa Array && V isa Array
            all(>(0), d) || throw(ArgumentError("marginal scales d must all be positive"))
            all(>(0), s) || throw(ArgumentError("direction scales s must all be positive"))
            # Zero-padded columns (used to hold the device rank slot open) carry s = 1
            # and contribute nothing; only active columns must be orthonormal.
            act = findall(!=(1.0), s)
            isempty(act) || opnorm(V[:, act]' * V[:, act] - I) < 1.0e-6 ||
                throw(ArgumentError("V must have orthonormal columns"))
        end
        return new{TB, TV_, TS}(b, d, V, s)
    end
end

function Base.show(io::IO, p::LowRankPreconditioner)
    return print(
        io, "LowRankPreconditioner(n = $(length(p.b)), rank = $(length(p.s)), " *
            "scales = $(round.(p.s; digits = 2)))"
    )
end

_affine_logdet(p::LowRankPreconditioner) = sum(log, p.d) + sum(log, p.s; init = 0.0)

function _affine_fwd(p::LowRankPreconditioner, z::AbstractVector)
    isempty(p.s) && return p.b .+ p.d .* z
    return p.b .+ p.d .* (z .+ p.V * ((p.s .- 1) .* (p.V' * z)))
end

function _affine_inv(p::LowRankPreconditioner, x::AbstractVector)
    w = (x .- p.b) ./ p.d
    isempty(p.s) && return w
    # (I + V (S - I) V')⁻¹ = I + V (S⁻¹ - I) V' since V is orthonormal.
    return w .+ p.V * ((inv.(p.s) .- 1) .* (p.V' * w))
end

# Pull a gradient back from the sampled space to the base-flat space. With x = b + A z the
# densities are related by a constant Jacobian, so ∇_z = Aᵀ ∇_x and ∇_x = A⁻ᵀ ∇_z. Writing
# A = Diagonal(d) * M with M = I + V (S - I) V' symmetric gives A⁻ᵀ = D⁻¹ M⁻¹ — the
# low-rank part acts BEFORE the diagonal here, the reverse of `_affine_inv`.
function _affine_invT(p::LowRankPreconditioner, g::AbstractVector)
    isempty(p.s) && return g ./ p.d
    return (g .+ p.V * ((inv.(p.s) .- 1) .* (p.V' * g))) ./ p.d
end

# First `k` columns of the thin Q factor of an n×m matrix (k ≤ m), without ever
# materializing the full n×n orthogonal factor. `Matrix(qr(A).Q)` builds that dense
# n×n via `orgqr`; for n ≈ 66000 that is a 35 GB allocation. Applying the compact Q
# to a thin n×k identity yields the same columns at O(n·k) cost.
function _thinq(A::AbstractMatrix, k::Integer)
    kk = min(Int(k), size(A, 2))
    return qr(A).Q * Matrix{eltype(A)}(I, size(A, 1), kk)
end

# Host mirrors for device-buffered stages. The sampler's host-side bookkeeping
# (logging the constrained draw each chunk, inverting positions for transfer) pushes
# plain host arrays through the transform; when its buffers live on the device they
# must be brought back to host for those calls, while traced/tracing calls must see
# the device buffers untouched. Dispatch on the INPUT array: plain host arrays get a
# hostified transform, everything else passes through.
_hostify(x::AbstractArray) = x isa Array ? x : Array(x)
_hostify(p::LowRankPreconditioner) =
    LowRankPreconditioner(_hostify(p.b), _hostify(p.d), _hostify(p.V), _hostify(p.s))
_devicebuffers(p::LowRankPreconditioner) = !(p.b isa Array)
_ishostvec(z) = z isa Array || (z isa SubArray && parent(z) isa Array)
_pre_for(p, z) = (_ishostvec(z) && _devicebuffers(p)) ? _hostify(p) : p

# The flat-transform node: the affine acts in LATENT space, before the inner transform
# (the mirror image of PT's `PushforwardTransform`, whose map acts after). The inner is
# the posterior's whole `asflat` node, so this consumes exactly `dimension(inner)`
# latent coordinates.
struct PreconditionedFlat{P, I <: TV.AbstractTransform} <: TV.VectorTransform
    pre::P
    inner::I
end

TV.dimension(t::PreconditionedFlat) = TV.dimension(t.inner)
PT.is_scalar_transport(t::PreconditionedFlat) = PT.is_scalar_transport(t.inner)

function TV.transform_with(
        flag::TV.LogJacFlag, t::PreconditionedFlat, z::AbstractVector, index
    )
    n = TV.dimension(t.inner)
    pre = _pre_for(t.pre, z)
    w = _affine_fwd(pre, view(z, index:(index + n - 1)))
    x, ℓi, _ = TV.transform_with(flag, t.inner, w, firstindex(w))
    flag isa TV.NoLogJac && return x, ℓi, index + n
    return x, ℓi + _affine_logdet(pre), index + n
end

TV.inverse_eltype(t::PreconditionedFlat, ::Type{T}) where {T} =
    TV.inverse_eltype(t.inner, T)

function TV.inverse_at!(z::AbstractVector, index, t::PreconditionedFlat, x)
    n = TV.dimension(t.inner)
    index′ = TV.inverse_at!(z, index, t.inner, x)
    zv = view(z, index:(index + n - 1))
    copyto!(zv, _affine_inv(_pre_for(t.pre, z), zv))
    return index′
end

# Lets `maybe_transport`/`resolve_disk_transport` accept the preconditioner as a latent
# "space": the sampled space is the flat one with the affine composed in front, and the
# result stays a `TransformedVLBIPosterior`, so every sampler/AD extension method
# dispatching on that concrete type keeps working.
_pre_dim(space::LowRankPreconditioner) = length(space.b)

function PT.transport_to(post::VLBIPosterior, space::LowRankPreconditioner)
    t0 = asflat(post).transform
    node0 = PT.transport_node(t0)
    n = TV.dimension(node0)
    _pre_dim(space) == n || throw(
        DimensionMismatch(
            "preconditioner dimension $(_pre_dim(space)) does not match the posterior " *
                "latent dimension $n — it was fit to a different model"
        )
    )
    node = PreconditionedFlat(space, node0)
    td = PT.TransportedDistribution(node, getfield(t0, :start), nothing)
    return TransformedVLBIPosterior(post, td)
end

# Convert a host-fitted preconditioner to device buffers, padding the low-rank stage to
# `rank_cap` columns (zero columns with `s = 1` are exact no-ops). Device buffers trace as
# runtime inputs of a compiled program, so a later refit updates them in place with
# `_update_device_pre!` — no recompile — as long as the shapes stay fixed. Host arrays
# instead enter a compiled program as constants, so refitting them forces one.
# Implemented by the Reactant extension.
function _device_pre(pre; kwargs...)
    return throw(
        ArgumentError("device preconditioner buffers need Reactant loaded; `using Reactant` first")
    )
end

# Copy a fresh host fit `h` into the live device buffers `dev`, whose shapes were fixed
# at `_device_pre` time. Implemented by the Reactant extension.
function _update_device_pre!(dev, h)
    return throw(
        ArgumentError("device preconditioner buffers need Reactant loaded; `using Reactant` first")
    )
end

# Caches reused across warmup refits: `draws`/`scores` map a warmup-store draw index to
# its base-flat position and flat-space score (both invariant across refits — the flat
# space never changes), and `fn` holds the score function built on first use. The windowed
# run seeds `draws` with the sampler's real base-flat draws as warmup produces them, so a
# refit reuses them and computes each score once; without the cache every refit rebuilds
# the score function and reloads the full window.
_new_refit_cache() = (;
    draws = Dict{Int, Vector{Float64}}(),
    scores = Dict{Int, Vector{Float64}}(),
    fn = Ref{Any}(nothing),
)

# Flat-space score of `post` evaluated on the device: `post` is moved to the Reactant
# device, and the returned closure maps a host `Vector{Float64}` to the host gradient of
# the flat log-density. The program is compiled once, at the shape of `x0`. Implemented by
# the Reactant + Enzyme extension.
function _compiled_flat_score(post, x0)
    return throw(
        ArgumentError(
            "a device flat score needs both Reactant and Enzyme loaded; `using Reactant, " *
                "Enzyme` first, or ask for the host score instead"
        )
    )
end

# Host score, through the posterior's own AD mode.
function _host_flat_score(post::VLBIPosterior)
    tflat = asflat(post)
    return x -> last(LogDensityProblems.logdensity_and_gradient(tflat, x))
end

# One flat-space score, through the cache's score function (built on first use).
function _flat_score!(cache, post::VLBIPosterior, x::Vector{Float64}; reactant::Bool)
    if cache.fn[] === nothing
        cache.fn[] = reactant ? _compiled_flat_score(post, x) : _host_flat_score(post)
    end
    return cache.fn[](x)
end

# nutpie-style first-score initialization (Seyboldt, Carlson & Carpenter 2026, §3.1):
# before ANY draws exist, one score evaluation at the start point sets the diagonal scales
# `d = 1/√|α⁰|` (the one-sample estimate of the conditional widths — the mass matrix
# `diag(|α⁰|)`, which also makes the sampler scale-free in the parameterization). The
# start point becomes the center, so a run passing the result as its `transport_method`
# begins warmup pre-scaled instead of on a unit metric; the first windowed refit then
# replaces this with the full Fisher fit. The score comes through the refit cache, so the
# score function is shared with every later refit.
_score_init_pre(post::VLBIPosterior, θ0; reactant::Bool) =
    _score_init_pre(post, θ0, _new_refit_cache(); reactant)

function _score_init_pre(post::VLBIPosterior, θ0, cache; reactant::Bool)
    x0 = inverse(asflat(post), θ0)
    g = _flat_score!(cache, post, x0; reactant)
    d = clamp.(inv.(sqrt.(abs.(g) .+ 1.0e-8)), 1.0e-4, 1.0e4)
    return LowRankPreconditioner(x0, d, zeros(length(x0), 0), Float64[])
end

# Symmetric-positive-definite geometric mean: the Σ solving Σ A Σ = B
# (Seyboldt, Carlson & Carpenter 2026, Algorithm 2).
function _spdm(A::Symmetric, B::Symmetric)
    FA = eigen(A)
    Ah = FA.vectors * Diagonal(sqrt.(FA.values)) * FA.vectors'
    Aih = FA.vectors * Diagonal(inv.(sqrt.(FA.values))) * FA.vectors'
    FM = eigen(Symmetric(Ah * B * Ah))
    Mh = FM.vectors * Diagonal(sqrt.(max.(FM.values, 0.0))) * FM.vectors'
    return Symmetric(Aih * Mh * Aih)
end

# The Fisher-divergence estimator of Seyboldt, Carlson & Carpenter (arXiv:2603.18845,
# Algorithm 1): fit the low-rank-plus-diagonal affine transform minimizing the sample
# Fisher divergence from the transformed density to a standard normal, using the flat
# draws `Z` AND their scores `G = ∇ log p` (n × N, same columns). Scores let the estimator
# bypass the Cramér–Rao limit of draw-only estimation — on the subspace spanned by draws
# and scores a Gaussian target is matched exactly — so short windows carry real
# information and no sampling-noise thresholds are needed there.
#
# The construction: diagonal scales are the coordinatewise geometric mean
# `σ = √(var(x)/var(α))` (marginal vs conditional width — the closed-form diagonal
# minimizer); the low-rank part comes from eigendecomposing the SPD geometric mean of the
# projected draw covariance and inverse score covariance on the joint draw-plus-score
# subspace, keeping eigenvalues `λ ≥ cutoff` or `≤ 1/cutoff` (wide AND stiff directions in
# one two-sided criterion, coupling included), up to `rank` by `|log λ|`. Returns a
# `LowRankPreconditioner` with `s = √λ`.
function _fisher_lowrank(
        Z::AbstractMatrix, G::AbstractMatrix;
        rank::Int, cutoff::Real = 2.0, γ::Real = 1.0e-5, s_floor::Real = 0.02, carry = nothing
    )
    n, N = size(Z)
    size(G) == (n, N) || throw(DimensionMismatch("draws and scores must align"))
    N >= 4 || throw(ArgumentError("need at least 4 draws, got $N"))
    vz = vec(var(Z; dims = 2))
    vg = vec(var(G; dims = 2))
    (all(>(0), vz) && all(>(0), vg)) || throw(
        ArgumentError("zero-variance coordinates in draws or scores")
    )
    σ = (vz ./ vg) .^ (1 // 4)                      # σ*² = √(var(x)/var(α))
    x̄ = vec(mean(Z; dims = 2))
    # The estimation is centered at the draw mean (Algorithm 1); the returned
    # transform's SHIFT is the score-informed center of Thm 2.2, μ* = x̄ + σ*² ⊙ ᾱ —
    # the mean score pulls it toward the mode, past the Cramér–Rao limit of x̄ alone.
    b = x̄ .+ σ .^ 2 .* vec(mean(G; dims = 2))
    X = (Z .- x̄) ./ σ                               # standardized draws
    A = (G .- mean(G; dims = 2)) .* σ               # standardized scores (contravariant)
    X0h = copy(X); A0h = copy(A)                    # pre-deflation, for carried re-scale
    # Accumulation (carry): carry the previous transform's directions forward so a
    # direction once found is never lost when a later window fails to see it. The carried
    # directions are re-expressed in the CURRENT standardization, deflated out of the
    # draws/scores so the residual fit spends its budget on NEW structure, and re-scaled
    # on the current draws — keeping the wider of carried and re-measured for wide
    # directions (a wide mode only ever gets wider as the chain frees).
    U0 = zeros(n, 0); s0 = Float64[]
    if carry !== nothing && !isempty(carry.s)
        U0 = _thinq((carry.d .* carry.V) ./ σ, length(carry.s))
        s0now = (vec(var(U0' * X0h; dims = 2)) ./ vec(var(U0' * A0h; dims = 2))) .^ (1 // 4)
        s0 = map((sc, sn) -> sc > 1 ? max(sc, sn) : sn, carry.s, s0now)
        X .-= U0 * (U0' * X)
        A .-= U0 * (U0' * A)
    end
    # Directions AND scales from ALL draws, using the in-sample geometric-mean
    # eigenvalues directly — matching nuts-rs (pymc-devs/nuts-rs,
    # transform/adapt/low_rank.rs): no cross-validation, no eigenvalue shrinkage. The γ
    # ridge on each covariance and the two-sided cutoff are the only regularization. The
    # mass-matrix eigenvalue λ is the scale²: s = √λ, wide for λ > cutoff, stiff for
    # λ < 1/cutoff.
    Q = _thinq(hcat(X, A), min(2 * (N - 1), n))
    Px = Q' * X
    Pa = Q' * A
    # cov = sample covariance + γI. nuts-rs uses XXᵀ/γ + I (an effective ridge ~γ/N); at
    # N ≪ n that weak ridge lets sampling noise in the joint draw+score subspace pass the
    # two-sided filter as spurious wide/stiff directions (dozens on ~150 draws). The
    # stronger γI directly suppresses them — the fitted low rank then tracks real
    # anisotropy, not noise. γ is a nuts-rs setting, so this is a different regularization
    # strength, not a different method.
    Cx = Symmetric(Px * Px' ./ (N - 1) + γ * I)
    Ca = Symmetric(Pa * Pa' ./ (N - 1) + γ * I)
    Σ = _spdm(Ca, Cx)                               # cov_draws # cov_grads⁻¹
    F = eigen(Σ)
    λ = F.values
    cand = findall(l -> l >= cutoff || l <= inv(cutoff), λ)
    if isempty(cand)
        isempty(U0) && return LowRankPreconditioner(b, σ, zeros(n, 0), Float64[])
        return LowRankPreconditioner(b, σ, hcat(U0), vcat(s0))
    end
    U = Q * F.vectors[:, cand]
    U = _thinq(U, length(cand))                     # re-orthonormalize after projection
    # In-sample geometric-mean scale s = √λ, floored on the stiff side. At N ≪ n a
    # direction's fitted scale can round toward zero (the geometric-mean eigenvalue
    # underflows); unfloored it compresses that axis by ~1/s and forces the step size
    # toward zero. The floor bounds that; real stiff directions sit well above it.
    s = clamp.(sqrt.(max.(λ[cand], eps())), s_floor, Inf)
    # `rank` bounds per-fit cost. A driver that pads device buffers to a fixed rank can
    # grow its cap to match, so this only bites at the configured ceiling.
    keep = sortperm(abs.(log.(s)); rev = true)
    keep = keep[1:min(max(rank - length(s0), 0), length(keep))]
    isempty(U0) && return LowRankPreconditioner(b, σ, U[:, keep], s[keep])
    # Deflating the carried span out of the draws/scores leaves U0 and the new
    # directions mutually orthogonal only up to roundoff, which trips the constructor's
    # orthonormality check. A column-order-preserving thin QR cleans the cross terms
    # without disturbing the per-direction scale assignment: the already-orthonormal
    # leading block reappears up to a per-column sign, and the scale is sign-invariant.
    V = _thinq(hcat(U0, U[:, keep]), length(s0) + length(keep))
    return LowRankPreconditioner(b, σ, V, vcat(s0, s[keep]))
end

# The preconditioner a pilot run itself sampled through, to be carried into a new fit as
# seed directions.
function _pilot_carry(root::AbstractString)
    tf = joinpath(root, "transport.jls")
    old = isfile(tf) ? deserialize(tf) : nothing
    old isa LowRankPreconditioner && return old
    @warn "augment requested but the pilot has no LowRankPreconditioner at $tf; " *
        "fitting without carried directions"
    return nothing
end

"""
    fit_preconditioner(pilot, post; rank=16, nsamples=2000, discard=0.0, augment=false) -> LowRankPreconditioner

Fit a [`LowRankPreconditioner`](@ref) for `post` from a pilot run's log by the
Fisher-divergence estimator of Seyboldt, Carlson & Carpenter (arXiv:2603.18845): draws
AND scores in flat space, one joint fit. `pilot` is an MCMC [`DiskStore`](@ref)
directory (a run's `outbase`, or `<outbase>/warmup` for the warmup log). The pilot must
have sampled the SAME model.

A run whose sampler adapted under [`FisherLowRank`](@ref) leaves a
`metric_adaptation.jls` checkpoint holding its warmup draws and scores already in
base-flat coordinates; that is used when present. Otherwise draws are reconstructed as
`inverse(asflat(post), θ)` from the stored chain and scored here — for a model with
angle-embedded `(sin, cos)` parameters that reconstruction returns the unit-radius
representative, collapsing each pair's radial direction into an artificial
zero-variance axis the chain never actually sampled.

Up to `nsamples` draws are used, thinned evenly after dropping the first `discard`
fraction of the chain; `rank` caps the number of corrected directions.

`augment = true` carries the pilot's own preconditioner (its `transport.jls`) into the
fit as seed directions, so successive rounds accumulate corrections instead of each
short log re-detecting only the widest and forgetting the rest.

`grad_reactant = true` scores the draws on the Reactant device, which requires both
Reactant and Enzyme to be loaded, instead of on the host.
"""
function fit_preconditioner(
        pilot::AbstractString, post::VLBIPosterior;
        rank::Int = 16, nsamples::Int = 2000, discard::Real = 0.0,
        augment::Bool = false, grad_reactant::Bool = false, refit_cache = nothing
    )
    0 <= discard < 1 || throw(ArgumentError("discard must be in [0, 1), got $discard"))
    # A run's own checkpoints live at its root; `pilot` may point at the warmup log one
    # level down.
    root = basename(abspath(pilot)) == "warmup" ? dirname(abspath(pilot)) : abspath(pilot)
    # The sampler's own base-flat draws and scores, when the run recorded them.
    af = joinpath(root, "metric_adaptation.jls")
    if isfile(af)
        st = deserialize(af)::FisherAdaptation
        N = length(st.draws)
        start = round(Int, discard * N) + 1
        sel = unique(round.(Int, range(start, N, length = min(192, N - start + 1))))
        pre = _fisher_lowrank(
            reduce(hcat, st.draws[sel]), reduce(hcat, st.scores[sel]);
            rank, carry = augment ? _pilot_carry(root) : nothing
        )
        @info "Fitted Fisher preconditioner from $(length(sel)) recorded draw/score pairs: $pre"
        return pre
    end
    carry = augment ? _pilot_carry(root) : nothing
    ntot = deserialize(joinpath(pilot, "parameters.jls")).params.nsamples
    start = round(Int, discard * ntot) + 1
    step = max(1, (ntot - start + 1) ÷ nsamples)
    usedidx = collect(start:step:ntot)
    tpost = asflat(post)
    n = dimension(tpost)
    # Fisher-divergence fit (Seyboldt+ 2026): draws + scores in flat space. Only the ≤192
    # columns used are loaded, each scored ONCE via the cache, and the score function is
    # built once — so a repeat fit against the same cache costs seconds.
    cache = isnothing(refit_cache) ? _new_refit_cache() : refit_cache
    sel = unique(round.(Int, range(1, length(usedidx), length = min(192, length(usedidx)))))
    need = usedidx[sel]
    Zf = zeros(n, length(need))
    Gf = zeros(n, length(need))
    for (k, idx) in enumerate(need)
        Zf[:, k] = get!(cache.draws, idx) do
            inverse(tpost, only(postsamples(load_samples(pilot, idx:idx))))
        end
        Gf[:, k] = get!(cache.scores, idx) do
            _flat_score!(cache, post, Zf[:, k]; reactant = grad_reactant)
        end
    end
    pre = _fisher_lowrank(Zf, Gf; rank, carry)
    @info "Fitted Fisher preconditioner from $(length(need)) draw/score pairs: $pre"
    return pre
end
