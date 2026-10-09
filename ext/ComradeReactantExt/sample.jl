using AbstractMCMC
using Serialization
using Random
using Printf
using Reactant: ProbProg

# ===========================================================================
# Sample-retention backends (reuse Comrade's MemoryStore / DiskStore configs)
# Both deal in transformed PosteriorSamples. A "sink" is driven by:
#   _open_sink -> create
#   _write_sink!(sink, store, chain, stats, state) -> persist one chunk
#   _close_sink -> finalize
# The sink takes only the data it needs to persist (chain, per-draw stats, state); the
# richer callback `info` is built separately in `sample_chunked` (mirroring the AdvancedHMC
# path, where serialization and the callback `info` are distinct steps).
# ===========================================================================

# --- MemoryStore: accumulate transformed chains + stats, build one PS at close ---
mutable struct _MemorySink
    chains::Vector{Any}
    stats::Vector{Any}
end
_open_sink(::MemoryStore, _tpost, _nsamples, _nscans, _stride; append::Bool = false) =
    _MemorySink(Any[], Any[])
function _write_sink!(sink::_MemorySink, ::MemoryStore, chain, stats, state)
    push!(sink.chains, chain)
    push!(sink.stats, stats)
    return nothing
end
function _close_sink(sink::_MemorySink, ::MemoryStore, metadata)
    chain = reduce(vcat, sink.chains)
    stats = map((cols...) -> reduce(vcat, cols), sink.stats...)
    return PosteriorSamples(chain, stats; metadata)
end

# --- DiskStore: write Comrade's `sample_to_disk` layout, return a DiskOutput ---
# On append (resume), `iter`/`nsamples_done` start from the existing chain so new
# chunks are numbered *after* the ones already on disk (matching AHMC's restart).
mutable struct _DiskSink
    outdir::String
    outbase::String
    stride::Int
    iter::Int
    nsamples_done::Int
end

# Scan files of a chain stored in `dir` in Comrade's on-disk layout (empty if none).
function _scan_files(dir::String)
    sampdir = joinpath(dir, "samples")
    isdir(sampdir) || return String[]
    return filter(f -> startswith(f, "output_scan_") && endswith(f, ".jls"), readdir(sampdir))
end

# Remove any chain already stored in `dir` (scan files, index, warmup log). A fresh
# (non-append) open must never write around stale scan files: the index is only rewritten
# when the first chunk lands, and until then `load_samples(dir)` would reassemble the old
# run's files into the new chain. Deleting them is reported with a warning.
function _clear_stale_chain(dir::String)
    scans = _scan_files(dir)
    pf = joinpath(dir, "parameters.jls")
    haveidx = isfile(pf)
    isempty(scans) && !haveidx && return nothing
    @warn "\"$dir\" already contains a sampled chain ($(length(scans)) scan file(s)$(haveidx ? " and an index" : "")); removing it to start fresh. Pass `restart = true` to continue an interrupted run instead."
    foreach(f -> rm(joinpath(dir, "samples", f)), scans)
    haveidx && rm(pf)
    return nothing
end

# Open Comrade's on-disk sample layout in `dir` (creating `dir/samples`). This is the one
# place the layout's open semantics live — the main chain and the warmup log both go
# through it. On append, new chunks are numbered after the ones already indexed; a fresh
# open clears any chain already there (`_clear_stale_chain`).
function _open_disk_layout(dir::String, stride::Int; append::Bool = false)
    sampdir = joinpath(dir, "samples")
    mkpath(sampdir)
    iter0, nsamples0 = 0, 0
    pf = joinpath(dir, "parameters.jls")
    if append && isfile(pf)
        prev = deserialize(pf).params
        iter0, nsamples0 = prev.nfiles, prev.nsamples
    elseif !append
        _clear_stale_chain(dir)
    end
    return _DiskSink(dir, joinpath(sampdir, "output_scan_"), stride, iter0, nsamples0)
end

_open_sink(store::DiskStore, _tpost, _nsamples, _nscans, stride; append::Bool = false) =
    _open_disk_layout(store.name, store.stride; append)

# Serialize one chunk of constrained draws + their column-oriented stats as the next scan
# file. The index is rewritten separately (`_write_index!`) so callers control the
# crash-consistency ordering around it.
function _write_scan_file!(sink::_DiskSink, chain, stats)
    sink.iter += 1
    sink.nsamples_done += length(chain)
    ps = PosteriorSamples(chain, stats)
    serialize(
        sink.outbase * Printf.@sprintf("%05d.jls", sink.iter),
        (samples = Comrade.postsamples(ps), stats = Comrade.samplerstats(ps))
    )
    return nothing
end

function _write_index!(sink::_DiskSink)
    out = Comrade.DiskOutput(abspath(sink.outdir), sink.iter, sink.stride, sink.nsamples_done)
    serialize(joinpath(sink.outdir, "parameters.jls"), (; params = out))
    return nothing
end

function _write_sink!(sink::_DiskSink, ::DiskStore, chain, stats, state)
    _write_scan_file!(sink, chain, stats)
    # resumable MCMC state checkpoint (latest wins)
    ProbProg.save_state(joinpath(sink.outdir, "state.jls"), state)
    # update parameters.jls every chunk so a crash mid-run leaves it current and
    # consistent with the files on disk (AHMC checkpoints its counter every batch).
    # The state checkpoint deliberately lands *before* the index: a crash between the two
    # leaves the index one chunk behind the state, and the restart then renumbers from the
    # index and overwrites the orphaned scan — never duplicating a segment.
    _write_index!(sink)
    return nothing
end
function _close_sink(sink::_DiskSink, store::DiskStore, metadata)
    out = Comrade.DiskOutput(abspath(store.name), sink.iter, sink.stride, sink.nsamples_done)
    serialize(joinpath(store.name, "parameters.jls"), (; params = out))
    # Persist the same metadata the MemoryStore path attaches to PosteriorSamples
    # (sampler tag, warmup_history, final_state, user metadata). Reload via load_samples.
    serialize(joinpath(store.name, "metadata.jls"), metadata)
    return out
end

# --- Warmup draw log: the adaptation chain, written as it goes ---
# ProbProg collects no per-step warmup trace (its warmup kernel runs with `num_samples = 0`,
# and asking for samples mid-warmup would inject non-adapting steps and break the
# bit-identical-to-fused-warmup property). What *is* exact and free is the state at the end
# of each warmup chunk, so that is what gets logged: one draw per chunk.
#
# The layout is Comrade's standard on-disk sample layout (`<dir>/samples/output_scan_%05d.jls`
# + `parameters.jls`) with `stride = 1`, so the warmup chain reads back with a plain
# `Comrade.load_samples(dir)`. The log IS a `_DiskSink` — same layout, same open/write
# helpers — so the two can never drift apart; a fresh open clears a warmup chain already
# in the directory (`_clear_stale_chain` via `_open_disk_layout`). The MCMC
# state is checkpointed separately (`<name>/state.jls` in `warmup_chunked`), so no state
# lands here.
_open_warmup_log(dir::String; append::Bool = false) = _open_disk_layout(dir, 1; append)

# `params` is one constrained draw; `stats` is a column-oriented NamedTuple of length-1
# vectors, matching the convention `_write_sink!` uses for the sampling chunks.
function _write_warmup_draw!(log::_DiskSink, params, stats)
    _write_scan_file!(log, [params], stats)
    # Rewrite the index every draw so an interrupted warmup still leaves a loadable chain.
    _write_index!(log)
    return nothing
end

"""
    _existing_disk_samples(store::DiskStore) -> (nfiles, nsamples)

How much of the chain is already on disk in `store.name` (0, 0 if none).
"""
function _existing_disk_samples(store::DiskStore)
    pf = joinpath(store.name, "parameters.jls")
    isfile(pf) || return (0, 0)
    prev = deserialize(pf).params
    return (prev.nfiles, prev.nsamples)
end

# ===========================================================================
# Host-side state view shared by both callback paths
# ===========================================================================

"""
    _current_state(state, tpost) -> NamedTuple

Plot-friendly, host-side view of a ProbProg `MCMCState`. Materializes the Reactant
arrays back to the host and transforms the current (unconstrained, flattened)
`position` into the constrained model parameters via `tpost`. This lets a callback
inspect or plot the *current* draw directly, e.g.

```julia
DiskStore(; name = "Results", callback = info -> plot(skymodel(post, info.params.sky)))
```

Fields:
  - `position::Vector`            current position in unconstrained, flattened space
  - `params`                      `transform(tpost, position)` — the constrained
                                  parameter `NamedTuple` (`(; sky[, instrument])`)
  - `potential_energy::Real`      potential energy at `position`
  - `gradient::Vector`            gradient of the log-density at `position`
  - `step_size::Real`             current leapfrog step size
  - `inverse_mass_matrix::Vector` current (diagonal) inverse mass matrix
"""
function _current_state(state, tpost)
    position = Array(state.position)
    return (;
        position,
        params = transform(tpost, vec(position)),
        potential_energy = only(Array(state.potential_energy)),
        gradient = Array(state.gradient),
        step_size = only(Array(state.step_size)),
        inverse_mass_matrix = Array(state.inverse_mass_matrix),
    )
end

# ===========================================================================
# Host-side moves between chunks
# ===========================================================================

# Run the `between_chunks` hook and hand back a state the next compiled kernel can take.
# A hook may return the position as a host array; it is moved to the device in the shape
# the kernels were compiled for. Whenever the position changed, the cached gradient and
# potential energy belong to the old point, so both are dropped and the next kernel
# recomputes them. The step size, metric, RNG and adaptation state carry over.
function _run_between_chunks(hook, state, tpost, info, host_rng)
    isnothing(hook) && return state
    before = Array(state.position)
    state = hook(state, tpost, info, host_rng)
    state isa ProbProg.MCMCState || throw(
        ArgumentError("`between_chunks` must return the MCMCState, got $(typeof(state))")
    )
    after = Array(state.position)
    length(after) == length(before) || throw(
        DimensionMismatch(
            "`between_chunks` changed the position length from $(length(before)) to $(length(after))"
        )
    )
    if after != before
        state.position = Reactant.to_rarray(reshape(convert(Array{eltype(before)}, after), size(before)))
        state.gradient = nothing
        state.potential_energy = nothing
    end
    return state
end

# `Comrade.MoveSet` on a device position. The base point and every map stay on the device;
# the preconditioner enters the maps as a runtime argument, so in-place refits of its buffers
# need no recompile. The log density is that of the base space over the device posterior, so
# no proposal goes through the preconditioner. The traceable moves' steps of a call run as one
# program (`_fused_steps`); other moves propose on the host.
function Comrade._move_kernels(ms::Comrade.MoveSet, tpost, position::Reactant.AbstractConcreteArray)
    zd = Reactant.to_rarray(collect(Float64, vec(Array(position))))
    pre() = Comrade._transport_pre(tpost)
    if isnothing(pre())
        fwd = finv = (p, z) -> z
    else
        fwd = Reactant.Compiler.compile((p, z) -> Comrade._affine_fwd(p, z), (pre(), zd))
        finv = Reactant.Compiler.compile((p, x) -> Comrade._affine_inv(p, x), (pre(), zd))
    end
    ctx = _device_move_context(ms, tpost)
    tbase = ctx.view.tbase
    ld = Reactant.Compiler.compile((tb, x) -> logdensityof(tb, x), (tbase, zd))
    # Distinct arguments: compiling with one array in both places traces them as one input.
    pd = Reactant.Compiler.compile(Comrade._prior_delta, (zd, copy(zd)))
    tj = Comrade._traced_moves(ms)
    fused = isempty(tj) ? nothing : _compile_fused(ms, tj, ctx, zd)
    function propose(j, x, step)
        x′, logdet = Comrade.propose(ms.moves[j], Array(x), step, ms.ctx)
        return Reactant.to_rarray(x′), Float64(logdet)
    end
    function run_fused(x, ℓ, steps, logu, active)
        R = size(steps, 1)
        x, ℓd, logα = fused(
            ctx, x, Reactant.ConcreteRNumber(ℓ), Reactant.to_rarray(vec(steps)),
            Reactant.to_rarray(vec(logu)), Reactant.to_rarray(vec(Float64.(active)))
        )
        return x, Float64(ℓd), reshape(Array(logα), R, :)
    end
    return (;
        load = p -> Reactant.to_rarray(collect(Float64, vec(Array(p)))),
        to_base = z -> fwd(pre(), z),
        from_base = x -> finv(pre(), x),
        logdensity = x -> Float64(ld(tbase, x)),
        prior_delta = (x, x′) -> Float64(pd(x, x′)),
        propose,
        fused = run_fused,
        download = Array,
        upload = (xh, x) -> Reactant.to_rarray(xh),
        store = (z, p) -> Reactant.to_rarray(reshape(convert(Array{eltype(p)}, Array(z)), size(p))),
    )
end

function _compile_fused(ms, tj, ctx, zd)
    moves = Tuple(ms.moves[tj])
    free = Tuple(ms.free[tj])
    R = maximum(ms.rounds[tj])
    v = Reactant.to_rarray(zeros(R * length(tj)))
    return Reactant.Compiler.compile(
        (c, x, ℓ, steps, logu, active) -> _fused_steps(c, moves, free, R, x, ℓ, steps, logu, active),
        (ctx, zd, Reactant.ConcreteRNumber(0.0), v, copy(v), copy(v))
    )
end

# The traceable moves' steps of one call as a traced loop over rounds, in the base space of
# `ctx.view`: step `(r, c)` (linear index `r + (c - 1) R`) is round `r` of `moves[c]`, run
# when `active > 0` and accepted when `logu < log α`; a NaN `log α` rejects. A move with
# `free[c]` is accepted without the likelihood (see `MoveSet`). Returns the final point,
# its log density and every step's `log α`.
function _fused_steps(ctx, moves, free, R, x, ℓ, steps, logu, active)
    logα = zero(steps)
    @trace track_numbers = false for r in 1:R
        x, ℓ, logα = _fused_round(ctx, moves, free, R, r, x, ℓ, steps, logu, active, logα)
    end
    return x, ℓ, logα
end

function _fused_round(ctx, moves, free, R, r, x, ℓ, steps, logu, active, logα)
    for c in eachindex(moves)
        i = r + (c - 1) * R
        x′, logdet = Comrade.propose(moves[c], x, Comrade._rget(steps, i), ctx)
        ℓ′ = free[c] ? ℓ + Comrade._prior_delta(x, x′) : logdensityof(ctx.view.tbase, x′)
        a = ℓ′ - ℓ + logdet
        accept = (Comrade._rget(active, i) > 0) & (Comrade._rget(logu, i) < a)
        x = ifelse.(accept, x′, x)
        ℓ = ifelse(accept, ℓ′, ℓ)
        Comrade.ComradeBase.rsetindex!(logα, a, i)
    end
    return x, ℓ, logα
end

# The move context over the device posterior of `tpost`, with the moves' context arrays on
# the device.
function _device_move_context(ms::Comrade.MoveSet, tpost)
    view = Comrade.CoordinateView(tpost.lpost, Comrade.space(ms.ctx.view))
    data = Base.structdiff(ms.ctx, NamedTuple{(:view,)})
    return merge((; view), map(v -> v isa AbstractArray ? Reactant.to_rarray(v) : v, data))
end

# ===========================================================================
# Default callbacks (called between rounds; return value is collected into history)
# ===========================================================================

"""
    default_warmup_callback(info) -> NamedTuple

Default `warmup_callback` on the `MemoryStore` path: log one line per warmup chunk and
return its progress, the current adapted step size, and the current draw (`params`).

Return values are collected into `warmup_history`, which ends up in `samplerinfo(out)`
(`MemoryStore`) or `metadata.jls` (`DiskStore`) — so keeping `params` here is what makes
the warmup chain inspectable on the `MemoryStore` path, where there is no disk to log to.
"""
function default_warmup_callback(info)
    @info "ReactantNUTS warmup" step = info.step total = info.total step_size = info.step_size
    return (; info.step, info.total, info.step_size, info.params)
end

"""
    default_warmup_callback_noparams(info) -> NamedTuple

Default `warmup_callback` on the `DiskStore` path: [`default_warmup_callback`](@ref)
without `params` in the return value. A `DiskStore` already streams every per-chunk draw
to `<name>/warmup` as it is produced (see [`warmup_chunked`](@ref)), so also retaining
each draw in `warmup_history` would hold the whole warmup chain in host memory for the
entire run and serialize it a second time into `metadata.jls` — for large models and a
fine `warmup_chunk`, potentially hundreds of MB of pure duplication.
"""
default_warmup_callback_noparams(info) =
    Base.structdiff(default_warmup_callback(info), NamedTuple{(:params,)})

# The post-warmup (sampling) callback is the shared `Comrade.default_disk_callback`, used
# for both the `MemoryStore` and `DiskStore` paths — see its docstring for the `info`
# fields. Warmup has its own ReactantNUTS-specific default (its `info` carries `num_warmup`).

# ===========================================================================
# Engine: single fused warmup + chunked sampling
# ===========================================================================

# The NUTS settings a compiled kernel needs. Kernels close over these values only, never the
# sampler: a closure is traced with everything it captures, and the metric adaptor may hold
# host data (e.g. a curvature function over a host posterior) that cannot be traced.
_nuts_settings(sampler::ReactantNUTS) = (;
    max_tree_depth = sampler.max_tree_depth,
    max_delta_energy = sampler.max_delta_energy,
    strong_zero = sampler.strong_zero,
)

# Compile a warmup kernel that advances an `MCMCState` by `nsteps` adaptation steps,
# threading the dual-averaging/Welford `adaptation` carried on the state. `total` and the
# runtime `warmup_offset` anchor Stan's windowed schedule to the *global* warmup length, so
# chopping warmup into chunks is bit-identical to one fused warmup (this is exactly what
# Reactant's own `run_chain` does — see EnzymeAD/Reactant.jl#2964). The offset is a runtime
# argument so one compiled kernel serves every chunk of a given length.
function _compile_warmup_kernel(
        state, ldf, tpost, nsteps::Int, total::Int, sampler::ReactantNUTS;
        adapt_mass_matrix::Bool = Comrade.adapts_welford(sampler.metric_adaptor)
    )
    nuts = _nuts_settings(sampler)
    fn = function (st::ProbProg.MCMCState, lf, off)
        # `_infer` returns (trace, diagnostics, log_densities, traced_result, state). The
        # `log_densities` slot was added in Reactant 0.2.275 — hence the compat lower bound.
        _, _, _, _, st_out = ProbProg._infer(
            st, lf, tpost;
            algorithm = :NUTS, num_warmup = nsteps, num_samples = 0,
            adapt_step_size = true, adapt_mass_matrix = adapt_mass_matrix,
            total_warmup = total, warmup_offset = off, nuts...,
        )
        return st_out
    end
    return Reactant.Compiler.compile(
        fn, (state, ldf, ConcreteRNumber(Int64(0))); optimize = :probprog
    )
end

"""
    warmup_chunked(rng, ldf, x0, tpost, sampler; chunk, callback, checkpoint,
                   progress_checkpoint, resume_state, warmup_done) -> (state, history)

Run Stan-windowed warmup over `sampler.n_adapts` steps, advancing the `MCMCState` (and the
dual-averaging/Welford `adaptation` it carries) in chunks of `chunk` steps. Each chunk runs
through a compiled `_infer` kernel with `total_warmup`/`warmup_offset` set to the *global*
warmup length and the steps already done, so the windowed schedule is unaffected by the
chunk boundaries — chunked warmup is bit-identical to one fused warmup.

`ldf(x, tpost)` is the log-density. Returns the post-warmup `MCMCState` and the vector of
per-chunk `callback` return values. The `info` passed to `callback` carries `step`, `total`,
`num_warmup` (== `total`), `step_size`, the host-side view (`position`/`params`/
`potential_energy`/`gradient`/`inverse_mass_matrix`, see [`_current_state`](@ref)), the
raw `state`, and `gradient_time`: the wall seconds of one gradient through the device
preconditioner, measured when a refit builds it (`nothing` before the first fit or
without Enzyme loaded).

If `checkpoint` is a path, the `MCMCState` is written there with `ProbProg.save_state` after
*every* chunk (it persists the adaptation accumulators), and the step count is recorded to
`progress_checkpoint`. Together these let `sample(...; restart=true)`
resume an interrupted warmup from the last completed chunk. To resume, pass the loaded state
as `resume_state` and the recorded step count as `warmup_done`.

`sampler.metric_adaptor` decides how the metric adapts. Under [`Comrade.FisherLowRank`](@ref)
this loop also drives the refits: at each scheduled step it fits a new
[`Comrade.LowRankPreconditioner`](@ref) from the base-flat draws and scores it has been
accumulating (both taken from the sampler's own state — no separate score evaluation and
no disk round-trip), swaps the latent space to it, and restarts dual averaging anchored at
the current step size. A refit updates the transform's device buffers in place where the
rank allows, so the compiled kernel is reused; a rank that outgrows the padded block
rebuilds them and costs one recompile. The accumulated draws are checkpointed to
`adaptation_checkpoint` alongside the state, and the fitted transform to
`transport_checkpoint`, so a restart resumes with the window it had built up.

`between_chunks(state, tpost, info, host_rng) -> state` runs at the end of every chunk,
after the refit, and may move `state.position` (e.g. a Metropolis–Hastings step); see the
`ReactantNUTS` `sample` method. Here `info` is `(; phase = :warmup, step, total, pre)`.

If `draws_dir` is a path, the draw at the end of each chunk is appended there as it is
produced, giving an inspectable warmup chain of `n_adapts ÷ chunk` draws:

```julia
warmup = Comrade.load_samples(joinpath("Results", "warmup"))   # PosteriorSamples
Comrade.samplerstats(warmup).step_size                         # adaptation trace
```

Per-draw stats are `step` (global warmup step), `step_size`, and `potential_energy`. The
cadence is the chunk length — ProbProg collects no per-step warmup trace, and forcing one
would perturb the chain (see the warmup-log note above `_open_warmup_log`), so use a smaller `chunk` for a
finer log. Passing `resume_state` appends to an existing log rather than restarting it.
"""
# The transform warmup started in, for an adaptor that keeps it (`carry = :keep`): host
# arrays, active columns only. Refits overwrite the transport checkpoint, so a fresh warmup
# also stores it beside that as `initial_transport.jls`, which a resumed warmup reads back.
function _initial_pre(adaptor, tpost, transport_checkpoint, fresh::Bool)
    (adaptor isa Comrade.FisherLowRank && adaptor.carry === :keep) || return nothing
    file = isnothing(transport_checkpoint) ? nothing :
        joinpath(dirname(transport_checkpoint), "initial_transport.jls")
    if fresh
        pre = Comrade._transport_pre(tpost)
        h = isnothing(pre) ? nothing : Comrade._active_directions(pre)
        isnothing(file) || serialize(file, h)
        return h
    end
    (isnothing(file) || !isfile(file)) && throw(
        ArgumentError("resuming a `carry = :keep` warmup needs its starting transform at $(something(file, "<no transport checkpoint>"))")
    )
    return deserialize(file)
end

function warmup_chunked(
        rng, ldf, x0, tpost, sampler::ReactantNUTS;
        chunk::Int, callback = default_warmup_callback,
        checkpoint = nothing, progress_checkpoint = nothing,
        adaptation_checkpoint = nothing, transport_checkpoint = nothing,
        draws_dir = nothing,
        resume_state = nothing, warmup_done::Int = 0,
        adaptation = nothing, segment_start::Int = 0,
        between_chunks = nothing, host_rng = Random.default_rng(),
    )

    na = sampler.n_adapts
    na > 0 || throw(ArgumentError("n_adapts must be positive"))
    chunk > 0 || throw(ArgumentError("warmup chunk length must be positive"))

    adaptor = sampler.metric_adaptor
    astate = isnothing(adaptation) ? Comrade.init_metric_adaptation(adaptor) : adaptation
    Comrade.check_metric_space(adaptor, Comrade._base_space(tpost))
    initial = _initial_pre(adaptor, tpost, transport_checkpoint, isnothing(resume_state))

    if isnothing(resume_state)
        T = eltype(x0)
        step0 = ConcreteRNumber(T(sampler.init_step_size))
        mass0 = Reactant.to_rarray(ones(T, length(x0)))
        # gradient/potential_energy start as `nothing` (computed on the first chunk) and there
        # are no adaptation accumulators yet. The `config` is left at its default — every NUTS
        # parameter is passed explicitly to `_infer`/`mcmc_logpdf`, so the stored config is
        # never read on this path.
        state = ProbProg.MCMCState(x0, nothing, nothing, step0, mass0, rng.seed, nothing)
        done = 0
    else
        state = resume_state
        done = warmup_done
    end

    # Append to the warmup chain when resuming, start a fresh one otherwise.
    wlog = isnothing(draws_dir) ? nothing :
        _open_warmup_log(draws_dir; append = !isnothing(resume_state))

    history = Any[]
    # Kernel cache: one compiled kernel per chunk length and state shape (gradient and
    # adaptation presence, position rank). A refit rebuilds the state with empty gradient
    # and adaptation slots and a between-chunks move empties the gradient slot, flipping
    # the shape back and forth — each shape must compile ONCE and be reused, not
    # recompiled per flip. A structural transform swap invalidates
    # the whole cache (the transform is baked into the kernels).
    kernels = Dict{Any, Any}()
    # Windowed-refit bookkeeping: `seg` anchors the Stan schedule to the CURRENT segment
    # (each refit restarts the schedule for the remaining warmup, since the geometry the
    # earlier windows adapted to no longer exists), and `pending` holds the refit steps
    # still ahead. `devpre` is the live padded device buffer block the transform reads
    # from, and `rankcap` its column count; a fit whose rank still fits is written into
    # it in place, so the compiled kernel stays valid.
    seg = segment_start
    pending = sort(filter(>(done), Comrade.metric_refit_steps(adaptor, na, chunk)))
    adapt_mm = Comrade.adapts_welford(adaptor) && seg == 0
    devpre = nothing
    rankcap = 0
    # One gradient's wall time through the device preconditioner, measured when it is built.
    gradient_time = nothing
    while done < na
        nsteps = min(chunk, na - done)
        !isempty(pending) && (nsteps = min(nsteps, first(pending) - done))
        # With the metric frozen (post-refit), the Stan window schedule is inert and
        # `total_warmup` need not track the segment — pinning it keeps the kernel cache
        # key stable across refits (the schedule anchor only matters while the Welford
        # metric adapts, i.e. in segment 0).
        total_c = adapt_mm ? (na - seg) : na
        off_c = adapt_mm ? (done - seg) : done
        key = (
            nsteps, isnothing(state.gradient), isnothing(state.adaptation),
            ndims(state.position), adapt_mm, total_c,
        )
        fresh = !haskey(kernels, key)
        kernel = get!(kernels, key) do
            _compile_warmup_kernel(
                state, ldf, tpost, nsteps, total_c, sampler;
                adapt_mass_matrix = adapt_mm
            )
        end
        # The kernel runs asynchronously; waiting for its result keeps the timing to the
        # NUTS steps alone.
        nuts_seconds = @elapsed begin
            state = kernel(state, ldf, ConcreteRNumber(Int64(off_c)))
            Reactant.synchronize(state.position)
        end
        done += nsteps

        cur = _current_state(state, tpost)
        # Record the draw and its score for the metric adaptor before anything is
        # checkpointed, so the checkpoint and the accumulated window agree on how much
        # warmup has been seen. `xbf` is this draw in base-flat coordinates, which is also
        # what a new transform has to be re-expressed from.
        xbf = Comrade.observe_draw!(
            adaptor, astate, Comrade._transport_pre(tpost), cur.position, cur.gradient
        )

        # Checkpoint AFTER the chunk so a crash resumes from completed work.
        if !isnothing(checkpoint)
            ProbProg.save_state(checkpoint, state)
            isnothing(adaptation_checkpoint) || isnothing(astate) ||
                serialize(adaptation_checkpoint, astate)
            isnothing(progress_checkpoint) ||
                serialize(
                progress_checkpoint,
                (; warmup_done = done, n_adapts = na, segment_start = seg)
            )
        end

        # Log the draw before the callback so a throwing callback can't lose it.
        isnothing(wlog) || _write_warmup_draw!(
            wlog, cur.params,
            (;
                step = [done],
                step_size = [cur.step_size],
                potential_energy = [cur.potential_energy],
            )
        )
        info = (;
            step = done, total = na, num_warmup = na,
            step_size = cur.step_size,
            position = cur.position, params = cur.params,
            potential_energy = cur.potential_energy, gradient = cur.gradient,
            inverse_mass_matrix = cur.inverse_mass_matrix,
            state, gradient_time, nsteps, nuts_seconds, fresh_kernel = fresh,
        )
        push!(history, callback(info))

        # Windowed refit. The adaptor fits a new latent space from the draws and scores
        # accumulated so far; the current draw is carried into it through base-flat, and
        # the sampler continues from there with dual averaging restarted at the current
        # adapted step size, the metric frozen at identity, and the windowed schedule
        # anchored to the new segment. Where the new rank still fits the padded device
        # block the buffers are overwritten in place and the compiled kernel is reused;
        # otherwise the block is rebuilt, which is a structural change and recompiles.
        if !isempty(pending) && done >= first(pending)
            popfirst!(pending)
            pre = Comrade.metric_refit(adaptor, astate; current = Comrade._transport_pre(tpost), initial, position = xbf)
            if isnothing(pre)
                @info "warmup metric refit at step $done skipped: too few draws recorded"
            else
                @info "warmup metric refit at step $done" pre
                base = Comrade._base_space(tpost)
                isnothing(transport_checkpoint) ||
                    serialize(transport_checkpoint, Comrade._in_space(base, pre))
                nrank = length(pre.s)
                if isnothing(devpre) || nrank > rankcap
                    grow = !isnothing(devpre)
                    rankcap = max(round(Int, 1.25 * nrank), 16)
                    devpre = Comrade._device_pre(pre; rank_cap = rankcap)
                    tpost = Comrade.maybe_transport(tpost.lpost, Comrade._in_space(base, devpre))
                    empty!(kernels)  # structural change: every cached kernel is stale
                    grow && @info "grew low-rank cap to $rankcap (fit rank $nrank); one recompile"
                    gradient_time = Comrade._gradient_seconds(tpost, Reactant.to_rarray(Comrade._affine_inv(pre, xbf)))
                    isnothing(gradient_time) ||
                        @info "one gradient through the device preconditioner: $(round(1000 * gradient_time; digits = 2)) ms"
                else
                    Comrade._update_device_pre!(devpre, pre)
                end
                T = eltype(cur.position)
                xnew = Reactant.to_rarray(T.(Comrade._affine_inv(pre, xbf)))
                step0 = ConcreteRNumber(T(cur.step_size))
                mass0 = Reactant.to_rarray(ones(T, length(cur.position)))
                state = ProbProg.MCMCState(
                    xnew, nothing, nothing, step0, mass0, state.rng, nothing
                )
                seg = done
                adapt_mm = false
                # No checkpoint here: `save_state` cannot serialize the rebuilt state's
                # empty gradient/adaptation slots. The next chunk's standard checkpoint
                # (one chunk later) persists the populated state together with the new
                # `segment_start`; until then the on-disk state predates the refit while
                # `transport_checkpoint` is already new, so a crash inside that short
                # window needs a fresh start rather than `restart = true`.
            end
        end

        state = _run_between_chunks(
            between_chunks, state, tpost,
            (; phase = :warmup, step = done, total = na, pre = Comrade._transport_pre(tpost)),
            host_rng,
        )
    end
    return state, history, tpost
end

"""
    sample_chunked(state, ldf, tpost, sampler::ReactantNUTS; num_samples, saveto, chunk_size)
        -> (state, out, history)

Draw `num_samples` post-warmup samples in chunks, threading `state` forward with
adaptation OFF (frozen metric + step size). Each chunk is transformed to constrained
space and routed to `saveto` (`MemoryStore` -> accumulate; `DiskStore` -> write
Comrade layout). `out` is a `PosteriorSamples` (memory) or `Comrade.DiskOutput`
(disk). After EVERY chunk the per-batch callback runs: `saveto.callback` for a `DiskStore`,
otherwise [`Comrade.default_disk_callback`](@ref). The `info` it receives is the same
`NamedTuple` documented under the `ReactantNUTS` `sample` method (common fields in
[`Comrade.default_disk_callback`](@ref), with the host-side `MCMCState` view from
[`_current_state`](@ref) under `info.extras`).

`between_chunks(state, tpost, info, host_rng) -> state` runs after the callback of every
chunk but the last, with `info = (; phase = :sampling, step, total, pre)`; see the
`ReactantNUTS` `sample` method.
"""
function sample_chunked(
        state, ldf, tpost, sampler::ReactantNUTS;
        num_samples::Int, saveto = MemoryStore(), chunk_size::Int = 100,
        append::Bool = false, metadata = Dict{Symbol, Any}(),
        between_chunks = nothing, host_rng = Random.default_rng(),
    )

    # The per-batch callback is configured solely on the DiskStore; MemoryStore just logs
    # with the default.
    callback = saveto isa DiskStore ? saveto.callback : Comrade.default_disk_callback

    num_samples > 0 || throw(ArgumentError("num_samples must be positive"))
    chunk = saveto isa DiskStore ? saveto.stride : chunk_size
    chunk > 0 || throw(ArgumentError("chunk length must be positive"))

    full, rem = divrem(num_samples, chunk)
    sizes = [fill(chunk, full); rem > 0 ? [rem] : Int[]]
    # Workaround for Reactant 0.2 mcmc_logpdf scalar-vs-1xi1 MLIR mismatch when
    # num_samples == 1: roll a trailing length-1 chunk into the previous one.
    if length(sizes) > 1 && sizes[end] == 1
        sizes[end - 1] += 1
        pop!(sizes)
    elseif length(sizes) == 1 && sizes[1] == 1
        throw(
            ArgumentError(
                "num_samples=1 currently hits a Reactant mcmc_logpdf MLIR bug; " *
                    "request at least 2 samples.",
            )
        )
    end
    nrounds = length(sizes)

    nuts = _nuts_settings(sampler)
    run_chunk = function (st, ns)
        return ProbProg.mcmc_logpdf(
            st, ldf, tpost;
            algorithm = :NUTS, num_warmup = 0, num_samples = ns,
            adapt_step_size = false, adapt_mass_matrix = false, nuts...
        )
    end

    # One compiled kernel per (chunk length, gradient presence): a between-chunks move
    # empties the gradient slot, which changes the kernel's inputs.
    compiled = Dict{Tuple{Int, Bool}, Any}()
    gradient_time = isnothing(Comrade._transport_pre(tpost)) ? nothing :
        Comrade._gradient_seconds(tpost, Reactant.to_rarray(vec(Array(state.position))))
    sink = _open_sink(saveto, tpost, num_samples, nrounds, chunk; append)
    history = Any[]
    ndone = 0

    for (round, ns) in enumerate(sizes)
        key = (ns, isnothing(state.gradient))
        fresh = !haskey(compiled, key)
        cfn = get!(compiled, key) do
            Reactant.Compiler.compile(run_chunk, (state, ns); optimize = :probprog)
        end
        t = @elapsed begin
            # (trace, diagnostics, log_densities, traced_result, state) — the
            # `log_densities` slot appeared in Reactant 0.2.275 (see the compat bound).
            samples, diagnostics, log_densities, _, state = cfn(state, ns)
            Reactant.synchronize(state.position)
        end

        raw = Array(samples)
        chain = [transform(tpost, r) for r in eachrow(raw)]
        # `diagnostics` is an `ns x 2` matrix (Reactant ≥ 0.2.275, see the compat bound).
        # Column 1 is a placeholder the impulse dialect reads `true` for every draw — even
        # for a chain frozen on a single point — so it carries no information. Column 2 is
        # the per-draw divergence signal: 0 for a healthy draw, 1 when the trajectory blows
        # past `max_delta_energy`, with intermediate values in between. The cross-backend
        # `numerical_error` contract (`default_disk_callback`) is a `Vector{Bool}` of
        # per-draw divergence flags — `count` is called on it — so any nonzero signal
        # marks the draw divergent.
        numerical_error = map(>(0), Array(diagnostics)[:, 2])

        # Persist the chunk first (serialization / accumulation), then build the callback
        # `info` and fire the callback — the same ordering the AdvancedHMC path uses.
        cur = _current_state(state, tpost)
        # Per-draw sampler statistics, one column each. The step size is frozen during
        # sampling and the wall time is measured per chunk, so both repeat over the chunk;
        # `time` divided by the cost of one gradient is the number of leapfrog steps per
        # draw, which ProbProg does not report. Compilation is outside `time`, but the first
        # run of a freshly compiled kernel (`extras.fresh_kernel`) pays one-time device
        # setup, so its `time` is an overestimate.
        stats = (;
            numerical_error,
            log_density = vec(Array(log_densities)),
            step_size = fill(cur.step_size, ns),
            time = fill(t / ns, ns),
        )
        _write_sink!(sink, saveto, chain, stats, state)

        info = (;
            round, nrounds, num_samples = ns, time = t,
            step_size = cur.step_size, params = cur.params, numerical_error,
            # Backend-specific, NOT part of the cross-backend contract — see `extras` in
            # `Comrade.default_disk_callback`. Here it's the host-side `MCMCState` view.
            extras = (;
                position = cur.position, gradient = cur.gradient,
                potential_energy = cur.potential_energy,
                inverse_mass_matrix = cur.inverse_mass_matrix,
                state, samples = raw, gradient_time, fresh_kernel = fresh,
            ),
        )
        push!(history, callback(info))
        ndone += ns

        if round < nrounds
            state = _run_between_chunks(
                between_chunks, state, tpost,
                (;
                    phase = :sampling, step = ndone, total = num_samples,
                    pre = Comrade._transport_pre(tpost),
                ),
                host_rng,
            )
        end
    end

    meta = merge(
        Dict{Symbol, Any}(:nsamples => num_samples, :sample_history => history),
        metadata,
    )
    return state, _close_sink(sink, saveto, meta), history
end

# ===========================================================================
# High-level entry point (AdvancedHMC-ext style)
# ===========================================================================

function _initial_position(rng, tpost, initial_params)
    x = if isnothing(initial_params)
        prior_sample(rng, tpost)
    else
        Comrade.inverse(tpost, initial_params)
    end
    return Reactant.to_rarray(x)
end

_default_ldf(x, tpost) = logdensityof(tpost, x)

"""
    sample(rng, post, sampler::ReactantNUTS, nsamples;
           transport_method=nothing,
           saveto=MemoryStore(), initial_params=nothing, restart=false,
           chunk_size=100, warmup_chunk=0, ldf=_default_ldf,
           host_rng=Random.default_rng(),
           warmup_callback=nothing, between_chunks=nothing)

Warm up (Stan-windowed adaptation run in chunks of the sampling size, or of `warmup_chunk`
if given — see [`warmup_chunked`](@ref)) then draw `nsamples` post-warmup samples from the Reactant
posterior `post`. Structured like the AdvancedHMC extension's `sample`, and
algorithmically identical to AdvancedHMC's NUTS adaptation.

`transport_method` picks the latent space (see `Comrade.maybe_transport`): the default
(`nothing`) flattens a raw `VLBIPosterior` with `asflat` and keeps an already-transformed
posterior as-is, while an explicit space (e.g. `TVFlat()`, `StdNormal()`) is always
honored, re-transporting if needed.

Returns a `NamedTuple` `(; out, state)` where `state` is the final ProbProg
`MCMCState` (held in memory, ready to inspect, plot via [`_current_state`](@ref), or
thread into a follow-up `sample_chunked`) and `out` is the standard Comrade output:

  - `saveto::MemoryStore` -> `out` is a `PosteriorSamples` (chain transformed to
    constrained space; `samplerstats` carries per-sample `numerical_error`,
    `log_density` (of the flat, unconstrained posterior), `step_size` and `time`
    (wall seconds per draw, averaged over the chunk);
    warmup/sample history + final state in the metadata).
  - `saveto::DiskStore` -> writes per-chunk `PosteriorSamples` to `saveto.name` in
    Comrade's on-disk layout; `out` is the `Comrade.DiskOutput` handle. Read the
    chain back with `Comrade.load_samples(out)` or `Comrade.load_samples(saveto.name)`.
    A resumable `MCMCState` is checkpointed to `<name>/state.jls` after every warmup
    chunk and after each sampling chunk (warmup progress in `<name>/warmup_progress.jls`).

## Metric adaptation

`sampler.metric_adaptor` selects how the metric adapts during warmup — see
[`Comrade.WelfordDiagonal`](@ref), [`Comrade.FixedMetric`](@ref) and
[`Comrade.FisherLowRank`](@ref). Under `FisherLowRank` the latent space is refit from the
run's own warmup draws and scores at scheduled steps and warmup continues in the fitted
coordinates; with a `DiskStore` each fit is written to `<name>/transport.jls` (which is
also what a restart resumes in) and the accumulated draws and scores to
`<name>/metric_adaptation.jls`, so a resumed warmup refits from the whole window rather
than only the steps taken since the restart.

## Inspecting warmup

The warmup (adaptation) chain is saved as it runs, so it can be looked at while the job is
still going or after a crash:

  - `saveto::DiskStore` -> one draw per warmup chunk is appended to `<name>/warmup` in the
    same on-disk layout as the main chain, so it loads with
    `Comrade.load_samples(joinpath(name, "warmup"))`. Per-draw `samplerstats` carry `step`
    (global warmup step), `step_size`, and `potential_energy` — e.g.
    `Comrade.samplerstats(warmup).step_size` is the step-size adaptation trace.
  - `saveto::MemoryStore` -> the same draws come back in `warmup_history` (each entry has
    `params`), reachable via `samplerinfo(out)[:warmup_history]`.

The log has one draw per warmup chunk, not per step: ProbProg does not collect a per-step
warmup trace, and making it do so would inject non-adapting steps into warmup and change
the chain. Set `warmup_chunk` for a finer log (it defaults to the sampling chunk size);
chunk length does not otherwise affect the chain.

## Callbacks

The per-batch callback is configured solely through the `DiskStore`: `saveto.callback` runs
after every batch when `saveto::DiskStore`, otherwise (for `MemoryStore`) the default
[`Comrade.default_disk_callback`](@ref) logger is used. The `info` it receives has the
common fields documented in [`Comrade.default_disk_callback`](@ref) plus an `extras`
`NamedTuple` of backend-specific data (not part of the cross-backend contract). On this
`ReactantNUTS` path `info.extras` is the host-side view of the current `MCMCState` (see
[`_current_state`](@ref)):

  - `extras.position`            : current draw in unconstrained, flattened space
  - `extras.gradient`            : gradient of the log-density at `position`
  - `extras.potential_energy`    : potential energy at `position`
  - `extras.inverse_mass_matrix` : current (diagonal) inverse mass matrix
  - `extras.state`               : the raw Reactant `MCMCState` (resumable checkpoint)
  - `extras.samples`             : the raw, unconstrained sample matrix for the batch
  - `extras.gradient_time`       : wall seconds of one gradient through the device
                                   preconditioner (`nothing` without one or without Enzyme)

`warmup_callback` runs once per warmup chunk; its `info` carries `step`/`total`/
`num_warmup` alongside the host-side state view (see [`default_warmup_callback`](@ref)
and [`warmup_chunked`](@ref)), plus `nsteps` (the NUTS steps of the chunk),
`nuts_seconds` (their wall time on the device, without compilation, moves or refits) and
`fresh_kernel` (whether the chunk ran a just-compiled kernel, whose first run pays one-time
device setup). The default (`nothing`) resolves per store: `MemoryStore`
keeps each chunk's draw in `warmup_history` ([`default_warmup_callback`](@ref)), while
`DiskStore` logs only scalars there ([`default_warmup_callback_noparams`](@ref)) since
its draws already stream to `<name>/warmup`. Overriding it replaces what lands in
`warmup_history`, so return `info.params` from a custom callback to keep the in-memory
warmup chain.

## Moves between chunks

`between_chunks(state, tpost, info, host_rng) -> state` runs on the host between compiled
NUTS chunks and may change `state.position`, e.g. with a Metropolis–Hastings step along a
direction NUTS explores slowly. Composing a posterior-preserving move with the NUTS chunks
preserves the posterior. During warmup it runs at the end of every warmup chunk (after the
callback and any metric refit); during sampling, after the callback of every chunk but
the last. `info` is `(; phase, step, total, pre)`: `phase` is `:warmup` or `:sampling`,
`step` counts the warmup steps (or post-warmup draws) done so far out of `total`, and `pre`
is the preconditioner composed into `tpost` (`nothing` for plain base-flat), so the
constrained point is `transform(tpost, vec(Array(state.position)))`. `host_rng` is the
`host_rng` keyword.

`state.position` is the latent point of `tpost`; the hook may replace it with a device or
host array of the same length. Whenever the position changes, its cached gradient and
potential energy are dropped and the next chunk recomputes them; the step size, metric,
and adaptation state carry over. Moves are not recorded: the chain holds the NUTS draws,
and the checkpointed `state.jls` is the state before the move, which a restart resumes
from. A hook keeps its own acceptance statistics. The warmup chunk length (by default the
`DiskStore` stride) sets how often it runs.

A fresh run (`restart=false`) requires a `DiskStore` directory with no chain in it:
sampling into a directory that already holds a previous run's chain or warmup log is an
error, so an old run is never silently overwritten.

`restart=true` (only with a `DiskStore`) continues an interrupted run, matching
AdvancedHMC's `restart`. Warmup is checkpointed after every chunk, so an interruption
*during* warmup resumes from the last completed chunk (the persisted adaptation state is
threaded back in); an interruption during sampling resumes from the last sampling chunk,
appending to the chain already on disk.

`nsamples` is the TOTAL target chain length, samples already on disk are counted,
and new chunks are numbered *after* them and appended (with `parameters.jls` grown
to the cumulative total). If the requested `nsamples` is already on disk, nothing is
drawn.
"""
function AbstractMCMC.sample(
        rng::Reactant.ReactantRNG, post, sampler::ReactantNUTS,
        nsamples::Int;
        transport_method = nothing,
        saveto = MemoryStore(), initial_params = nothing, restart::Bool = false,
        chunk_size::Int = 100, warmup_chunk::Int = 0,
        ldf = _default_ldf, host_rng = Random.default_rng(),
        warmup_callback = nothing, between_chunks = nothing
    )

    # Default warmup callback: `MemoryStore` keeps the per-chunk draw in `warmup_history`
    # (its only record of the warmup chain); `DiskStore` drops it there because the same
    # draws already stream to `<name>/warmup`, and duplicating them in host memory +
    # `metadata.jls` is pure waste. An explicit `warmup_callback` always wins.
    if isnothing(warmup_callback)
        warmup_callback = saveto isa DiskStore ?
            default_warmup_callback_noparams : default_warmup_callback
    end

    # Persist/reload the latent space (DiskStore only) so a restart resumes in the same space
    # it was launched in instead of silently defaulting; MemoryStore just honors the kwarg.
    tpost = _device_transport(
        Comrade.resolve_disk_transport(
            post, saveto isa DiskStore ? saveto.name : nothing, restart, transport_method
        )
    )

    # Checkpoint paths (DiskStore only): the resumable MCMCState and the warmup step counter.
    # Warmup is chunked at the same size as sampling.
    state_ckpt, progress_ckpt, adapt_ckpt, transport_ckpt = if saveto isa DiskStore
        mkpath(saveto.name)
        (
            joinpath(saveto.name, "state.jls"),
            joinpath(saveto.name, "warmup_progress.jls"),
            joinpath(saveto.name, "metric_adaptation.jls"),
            joinpath(saveto.name, "transport.jls"),
        )
    else
        (nothing, nothing, nothing, nothing)
    end
    # Warmup chunking defaults to the sampling chunk size. It also sets the cadence of the
    # warmup draw log (one draw per chunk), so `warmup_chunk` lets that be made finer
    # without touching the sampling stride. Chunk size does not affect the chain itself —
    # the windowed schedule is anchored to `total_warmup`/`warmup_offset`.
    warmup_chunk = if warmup_chunk > 0
        warmup_chunk
    else
        saveto isa DiskStore ? saveto.stride : chunk_size
    end
    # DiskStore: log the adaptation chain to `<name>/warmup` as it goes.
    warmup_draws_dir = saveto isa DiskStore ? joinpath(saveto.name, "warmup") : nothing

    if restart
        saveto isa DiskStore ||
            throw(ArgumentError("restart=true requires saveto::DiskStore"))
        isfile(state_ckpt) ||
            throw(ArgumentError("cannot restart: no state checkpoint at $state_ckpt"))
    elseif saveto isa DiskStore
        # A fresh run clears a previous run's chain and warmup log up front (see
        # `_clear_stale_chain`) so the warning lands before the multi-minute warmup compile.
        _clear_stale_chain(saveto.name)
        _clear_stale_chain(warmup_draws_dir)
    end

    if !restart
        x0 = _initial_position(host_rng, tpost, initial_params)
        @info "ReactantNUTS warmup" n_adapts = sampler.n_adapts chunk = warmup_chunk
        state, warmup_history, tpost = warmup_chunked(
            rng, ldf, x0, tpost, sampler;
            chunk = warmup_chunk, callback = warmup_callback,
            checkpoint = state_ckpt, progress_checkpoint = progress_ckpt,
            adaptation_checkpoint = adapt_ckpt, transport_checkpoint = transport_ckpt,
            draws_dir = warmup_draws_dir, between_chunks, host_rng,
        )
    else
        # Restart: load the checkpoint. If warmup did not finish (recorded step count <
        # n_adapts) resume it from the last completed chunk, threading the persisted
        # adaptation accumulators; otherwise skip straight to sampling. `warmup_progress.jls`
        # is absent for already-completed warmups, so default to "complete".
        state = ProbProg.load_state(state_ckpt)
        prog = isfile(progress_ckpt) ? deserialize(progress_ckpt) : nothing
        warmup_done = isnothing(prog) ? sampler.n_adapts : prog.warmup_done
        # Refit segments: resume inside the segment the checkpoint belongs to. The
        # transform itself was re-persisted by the refit callable, so `tpost` above is
        # already the segment's transform.
        seg0 = (isnothing(prog) || !haskey(pairs(prog), :segment_start)) ? 0 :
            prog.segment_start
        # The accumulated warmup draws and scores are the metric adaptor's window; without
        # them a resumed run would refit from only the steps taken after the restart.
        adaptation = isfile(adapt_ckpt) ? deserialize(adapt_ckpt) : nothing
        if warmup_done < sampler.n_adapts
            @info "ReactantNUTS restart: resuming warmup" done = warmup_done n_adapts = sampler.n_adapts
            state, warmup_history, tpost = warmup_chunked(
                rng, ldf, nothing, tpost, sampler;
                chunk = warmup_chunk, callback = warmup_callback,
                checkpoint = state_ckpt, progress_checkpoint = progress_ckpt,
                adaptation_checkpoint = adapt_ckpt, transport_checkpoint = transport_ckpt,
                draws_dir = warmup_draws_dir,
                resume_state = state, warmup_done = warmup_done,
                adaptation = adaptation, segment_start = seg0,
                between_chunks, host_rng,
            )
        else
            @info "ReactantNUTS restart: warmup complete, skipping" dir = saveto.name
            warmup_history = nothing
        end
    end

    # On restart, `nsamples` is the TOTAL target: append only what's left on disk.
    remaining = nsamples
    if restart
        _, done = _existing_disk_samples(saveto)
        remaining = nsamples - done
        if remaining <= 0
            @warn "Requested $nsamples samples but $done already on disk; nothing to do."
            out = deserialize(joinpath(saveto.name, "parameters.jls")).params
            return (; out, state)
        end
        @info "Appending $remaining samples to existing chain of $done" total = nsamples
    end

    metadata = Dict{Symbol, Any}(
        :sampler => :ReactantNUTS,
        :warmup_history => warmup_history,
        :final_state => state,
    )
    state, out, _ = sample_chunked(
        state, ldf, tpost, sampler;
        num_samples = remaining, saveto, chunk_size, append = restart, metadata,
        between_chunks, host_rng,
    )
    # sample_history is already merged into metadata by sample_chunked, so it lives
    # in samplerinfo(out) for MemoryStore and in metadata.jls for DiskStore.
    return (; out, state)
end

"""
    sample(post, sampler::ReactantNUTS, nsamples; kwargs...)

Convenience method that builds a default `ReactantRNG` and forwards to the
main `sample`.
"""
function AbstractMCMC.sample(post, sampler::ReactantNUTS, nsamples::Int; kwargs...)
    rng = Reactant.ReactantRNG()
    return AbstractMCMC.sample(rng, post, sampler, nsamples; kwargs...)
end
