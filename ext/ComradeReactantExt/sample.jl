using AbstractMCMC
using Serialization
using Random
using Printf
using Reactant: ProbProg

# ===========================================================================
# Sample-retention backends (reuse Comrade's MemoryStore / DiskStore configs)
# Both deal in transformed PosteriorSamples. A "sink" is driven by:
#   _open_sink -> create
#   _write_sink!(sink, store, chain, numerical_error, state) -> persist one chunk
#   _close_sink -> finalize
# The sink takes only the data it needs to persist (chain, numerical_error, state); the
# richer callback `info` is built separately in `sample_chunked` (mirroring the AdvancedHMC
# path, where serialization and the callback `info` are distinct steps).
# ===========================================================================

# --- MemoryStore: accumulate transformed chains + stats, build one PS at close ---
mutable struct _MemorySink
    chains::Vector{Any}
    nerr::Vector{Any}
end
_open_sink(::MemoryStore, _tpost, _nsamples, _nscans, _stride; append::Bool = false) =
    _MemorySink(Any[], Any[])
function _write_sink!(sink::_MemorySink, ::MemoryStore, chain, numerical_error, state)
    push!(sink.chains, chain)
    push!(sink.nerr, numerical_error)
    return nothing
end
function _close_sink(sink::_MemorySink, ::MemoryStore, metadata)
    chain = reduce(vcat, sink.chains)
    stats = (; numerical_error = reduce(vcat, sink.nerr))
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

# Open Comrade's on-disk sample layout in `dir` (creating `dir/samples`). This is the one
# place the layout's open semantics live — the main chain and the warmup log both go
# through it. On append, new chunks are numbered after the ones already indexed. On a
# fresh open, any chain already in `dir` is dropped first — the scan files (only ones
# matching our own naming) *and* the index. The index must go with the files: it is only
# rewritten when the first chunk lands (behind the multi-minute kernel compile for the
# warmup log), and until then `load_samples(dir)` — an advertised mid-run workflow —
# would read the old run's index and crash on (or splice in) the files deleted here. It
# also keeps a later `append` open from adopting the old run's counter as its start.
function _open_disk_layout(dir::String, stride::Int; append::Bool = false)
    sampdir = joinpath(dir, "samples")
    mkpath(sampdir)
    iter0, nsamples0 = 0, 0
    pf = joinpath(dir, "parameters.jls")
    if append && isfile(pf)
        prev = deserialize(pf).params
        iter0, nsamples0 = prev.nfiles, prev.nsamples
    elseif !append
        for f in readdir(sampdir)
            startswith(f, "output_scan_") && endswith(f, ".jls") &&
                rm(joinpath(sampdir, f))
        end
        rm(pf; force = true)
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

function _write_sink!(sink::_DiskSink, ::DiskStore, chain, numerical_error, state)
    _write_scan_file!(sink, chain, (; numerical_error))
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
# helpers — so the two can never drift apart; a fresh open drops any previous warmup
# chain (files + index) via `_open_disk_layout`. The MCMC state is checkpointed
# separately (`<name>/state.jls` in `warmup_chunked`), so no state lands here.
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

# Compile a warmup kernel that advances an `MCMCState` by `nsteps` adaptation steps,
# threading the dual-averaging/Welford `adaptation` carried on the state. `total` and the
# runtime `warmup_offset` anchor Stan's windowed schedule to the *global* warmup length, so
# chopping warmup into chunks is bit-identical to one fused warmup (this is exactly what
# Reactant's own `run_chain` does — see EnzymeAD/Reactant.jl#2964). The offset is a runtime
# argument so one compiled kernel serves every chunk of a given length.
function _compile_warmup_kernel(state, ldf, tpost, nsteps::Int, total::Int, sampler::ReactantNUTS)
    fn = function (st::ProbProg.MCMCState, lf, off)
        # `_infer` returns (trace, diagnostics, log_densities, traced_result, state). The
        # `log_densities` slot was added in Reactant 0.2.275 — hence the compat lower bound.
        _, _, _, _, st_out = ProbProg._infer(
            st, lf, tpost;
            algorithm = :NUTS, num_warmup = nsteps, num_samples = 0,
            adapt_step_size = true, adapt_mass_matrix = true,
            total_warmup = total, warmup_offset = off,
            max_tree_depth = sampler.max_tree_depth,
            max_delta_energy = sampler.max_delta_energy,
            strong_zero = sampler.strong_zero,
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
`potential_energy`/`gradient`/`inverse_mass_matrix`, see [`_current_state`](@ref)), and the
raw `state`.

If `checkpoint` is a path, the `MCMCState` is written there with `ProbProg.save_state` after
*every* chunk (it persists the adaptation accumulators), and the step count is recorded to
`progress_checkpoint`. Together these let `sample(...; restart=true)`
resume an interrupted warmup from the last completed chunk. To resume, pass the loaded state
as `resume_state` and the recorded step count as `warmup_done`.

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
function warmup_chunked(
        rng, ldf, x0, tpost, sampler::ReactantNUTS;
        chunk::Int, callback = default_warmup_callback,
        checkpoint = nothing, progress_checkpoint = nothing,
        draws_dir = nothing,
        resume_state = nothing, warmup_done::Int = 0,
    )

    na = sampler.n_adapts
    na > 0 || throw(ArgumentError("n_adapts must be positive"))
    chunk > 0 || throw(ArgumentError("warmup chunk length must be positive"))

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
    kernel = nothing
    kernel_key = nothing
    while done < na
        nsteps = min(chunk, na - done)
        # Recompile when the chunk length changes or when the state gains its gradient (and
        # adaptation) after the very first chunk — mirrors Reactant's `_run_chain_chunked`.
        key = (nsteps, isnothing(state.gradient))
        if key != kernel_key
            kernel = _compile_warmup_kernel(state, ldf, tpost, nsteps, na, sampler)
            kernel_key = key
        end
        state = kernel(state, ldf, ConcreteRNumber(Int64(done)))
        done += nsteps

        # Checkpoint AFTER the chunk so a crash resumes from completed work.
        if !isnothing(checkpoint)
            ProbProg.save_state(checkpoint, state)
            isnothing(progress_checkpoint) ||
                serialize(progress_checkpoint, (; warmup_done = done, n_adapts = na))
        end

        cur = _current_state(state, tpost)
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
            state,
        )
        push!(history, callback(info))
    end
    return state, history
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
`NamedTuple` documented under the `ReactantNUTS` [`sample`](@ref) method (common fields in
[`Comrade.default_disk_callback`](@ref), with the host-side `MCMCState` view from
[`_current_state`](@ref) under `info.extras`).
"""
function sample_chunked(
        state, ldf, tpost, sampler::ReactantNUTS;
        num_samples::Int, saveto = MemoryStore(), chunk_size::Int = 100,
        append::Bool = false, metadata = Dict{Symbol, Any}()
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

    run_chunk = function (st, ns)
        return ProbProg.mcmc_logpdf(
            st, ldf, tpost;
            algorithm = :NUTS, num_warmup = 0, num_samples = ns,
            max_tree_depth = sampler.max_tree_depth,
            max_delta_energy = sampler.max_delta_energy,
            adapt_step_size = false, adapt_mass_matrix = false,
            strong_zero = sampler.strong_zero
        )
    end

    compiled = Dict{Int, Any}()
    sink = _open_sink(saveto, tpost, num_samples, nrounds, chunk; append)
    history = Any[]

    for (round, ns) in enumerate(sizes)
        cfn = get!(compiled, ns) do
            Reactant.Compiler.compile(run_chunk, (state, ns); optimize = :probprog)
        end
        t = @elapsed begin
            # (trace, diagnostics, log_densities, traced_result, state) — the
            # `log_densities` slot appeared in Reactant 0.2.275 (see the compat bound).
            samples, diagnostics, _, _, state = cfn(state, ns)
        end

        raw = Array(samples)
        chain = [transform(tpost, r) for r in eachrow(raw)]
        # Reactant 0.2.275 widened `diagnostics` from a length-`ns` Bool vector to an `ns x 2`
        # matrix. Column 1 is the old vector: the impulse dialect calls it a "placeholder for
        # future expansion" and it reads `true` for every draw — including a chain frozen on a
        # single point — so `.!col1` (what this used to compute) was uniformly `false` and
        # never flagged anything. Column 2 is the actual per-draw divergence flag: it is 0 for
        # a healthy chain, 1 when the trajectory blows past `max_delta_energy`, and takes
        # intermediate values in between. That matches `numerical_error`'s cross-backend
        # meaning (AdvancedHMC reports divergences here), so it is used directly, unnegated.
        numerical_error = Array(diagnostics)[:, 2]

        # Persist the chunk first (serialization / accumulation), then build the callback
        # `info` and fire the callback — the same ordering the AdvancedHMC path uses.
        _write_sink!(sink, saveto, chain, numerical_error, state)

        cur = _current_state(state, tpost)
        info = (;
            round, nrounds, num_samples = ns, time = t,
            step_size = cur.step_size, params = cur.params, numerical_error,
            # Backend-specific, NOT part of the cross-backend contract — see `extras` in
            # `Comrade.default_disk_callback`. Here it's the host-side `MCMCState` view.
            extras = (;
                position = cur.position, gradient = cur.gradient,
                potential_energy = cur.potential_energy,
                inverse_mass_matrix = cur.inverse_mass_matrix,
                state, samples = raw,
            ),
        )
        push!(history, callback(info))
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
           warmup_callback=nothing)

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
    constrained space; `samplerstats` carries per-sample `numerical_error`;
    warmup/sample history + final state in the metadata).
  - `saveto::DiskStore` -> writes per-chunk `PosteriorSamples` to `saveto.name` in
    Comrade's on-disk layout; `out` is the `Comrade.DiskOutput` handle. Read the
    chain back with `Comrade.load_samples(out)` or `Comrade.load_samples(saveto.name)`.
    A resumable `MCMCState` is checkpointed to `<name>/state.jls` after every warmup
    chunk and after each sampling chunk (warmup progress in `<name>/warmup_progress.jls`).

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

`warmup_callback` runs once per warmup chunk; its `info` carries `step`/`total`/
`num_warmup` alongside the host-side state view (see [`default_warmup_callback`](@ref)
and [`warmup_chunked`](@ref)). The default (`nothing`) resolves per store: `MemoryStore`
keeps each chunk's draw in `warmup_history` ([`default_warmup_callback`](@ref)), while
`DiskStore` logs only scalars there ([`default_warmup_callback_noparams`](@ref)) since
its draws already stream to `<name>/warmup`. Overriding it replaces what lands in
`warmup_history`, so return `info.params` from a custom callback to keep the in-memory
warmup chain.

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
        warmup_callback = nothing
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
    tpost = Comrade.resolve_disk_transport(
        post, saveto isa DiskStore ? saveto.name : nothing, restart, transport_method
    )

    # Checkpoint paths (DiskStore only): the resumable MCMCState and the warmup step counter.
    # Warmup is chunked at the same size as sampling.
    state_ckpt, progress_ckpt = if saveto isa DiskStore
        mkpath(saveto.name)
        (joinpath(saveto.name, "state.jls"), joinpath(saveto.name, "warmup_progress.jls"))
    else
        (nothing, nothing)
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
    end

    if !restart
        x0 = _initial_position(host_rng, tpost, initial_params)
        @info "ReactantNUTS warmup" n_adapts = sampler.n_adapts chunk = warmup_chunk
        state, warmup_history = warmup_chunked(
            rng, ldf, x0, tpost, sampler;
            chunk = warmup_chunk, callback = warmup_callback,
            checkpoint = state_ckpt, progress_checkpoint = progress_ckpt,
            draws_dir = warmup_draws_dir,
        )
    else
        # Restart: load the checkpoint. If warmup did not finish (recorded step count <
        # n_adapts) resume it from the last completed chunk, threading the persisted
        # adaptation accumulators; otherwise skip straight to sampling. `warmup_progress.jls`
        # is absent for already-completed warmups, so default to "complete".
        state = ProbProg.load_state(state_ckpt)
        warmup_done = isfile(progress_ckpt) ?
            deserialize(progress_ckpt).warmup_done : sampler.n_adapts
        if warmup_done < sampler.n_adapts
            @info "ReactantNUTS restart: resuming warmup" done = warmup_done n_adapts = sampler.n_adapts
            state, warmup_history = warmup_chunked(
                rng, ldf, nothing, tpost, sampler;
                chunk = warmup_chunk, callback = warmup_callback,
                checkpoint = state_ckpt, progress_checkpoint = progress_ckpt,
                draws_dir = warmup_draws_dir,
                resume_state = state, warmup_done = warmup_done,
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
        num_samples = remaining, saveto, chunk_size, append = restart, metadata
    )
    # sample_history is already merged into metadata by sample_chunked, so it lives
    # in samplerinfo(out) for MemoryStore and in metadata.jls for DiskStore.
    return (; out, state)
end

"""
    sample(post, sampler::ReactantNUTS, nsamples; kwargs...)

Convenience method that builds a default `ReactantRNG` and forwards to the
main [`sample`](@ref).
"""
function AbstractMCMC.sample(post, sampler::ReactantNUTS, nsamples::Int; kwargs...)
    rng = Reactant.ReactantRNG()
    return AbstractMCMC.sample(rng, post, sampler, nsamples; kwargs...)
end
