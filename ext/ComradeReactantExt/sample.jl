using AbstractMCMC
using Serialization
using Random
using HypercubeTransform
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
function _open_sink(store::DiskStore, _tpost, _nsamples, _nscans, stride; append::Bool = false)
    mkpath(joinpath(store.name, "samples"))
    iter0, nsamples0 = 0, 0
    pf = joinpath(store.name, "parameters.jls")
    if append && isfile(pf)
        prev = deserialize(pf).params
        iter0, nsamples0 = prev.nfiles, prev.nsamples
    end
    return _DiskSink(
        store.name, joinpath(store.name, "samples", "output_scan_"),
        store.stride, iter0, nsamples0
    )
end
function _write_sink!(sink::_DiskSink, ::DiskStore, chain, numerical_error, state)
    sink.iter += 1
    sink.nsamples_done += length(chain)
    ps = PosteriorSamples(chain, (; numerical_error))
    serialize(
        sink.outbase * Printf.@sprintf("%05d.jls", sink.iter),
        (samples = Comrade.postsamples(ps), stats = Comrade.samplerstats(ps))
    )
    # resumable MCMC state checkpoint (latest wins)
    ProbProg.save_state(joinpath(sink.outdir, "state.jls"), state)
    # update parameters.jls every chunk so a crash mid-run leaves it current and
    # consistent with the files on disk (AHMC checkpoints its counter every batch).
    out = Comrade.DiskOutput(abspath(sink.outdir), sink.iter, sink.stride, sink.nsamples_done)
    serialize(joinpath(sink.outdir, "parameters.jls"), (; params = out))
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

Default `warmup_callback`: log one line after warmup and return the adapted step size.
"""
function default_warmup_callback(info)
    @info "ReactantNUTS warmup complete" num_warmup = info.num_warmup step_size = info.step_size
    return (; info.num_warmup, info.step_size)
end

# The post-warmup (sampling) callback is the shared `Comrade.default_disk_callback`, used
# for both the `MemoryStore` and `DiskStore` paths — see its docstring for the `info`
# fields. Warmup has its own ReactantNUTS-specific default (its `info` carries `num_warmup`).

# ===========================================================================
# Engine: single fused warmup + chunked sampling
# ===========================================================================

"""
    warmup_single(rng, ldf, x0, tpost, sampler::ReactantNUTS; callback, checkpoint)
        -> (state, history)

Run warmup as a **single** `mcmc_logpdf` call over `sampler.n_adapts` steps. The
ProbProg NUTS pass performs Stan's windowed adaptation internally (init/term buffers,
doubling windows, one continuous dual-averaging step size + Welford diagonal mass
matrix), which is algorithmically identical to AdvancedHMC's `StanHMCAdaptor`. This is
deliberately *not* re-windowed on the host: chopping warmup into per-window
`mcmc_logpdf` calls restarts dual averaging (and re-inflates its prox-center to
`10*ϵ`) every window, so the step size never converges and the chain diverges.

`ldf(x, tpost)` is the log-density. Returns the post-warmup `MCMCState` and the
single `callback` return value. The `info` NamedTuple passed to `callback` has fields:
  num_warmup, step_size, position, params, potential_energy, gradient,
  inverse_mass_matrix, state, samples, diagnostics.
`state` is the raw Reactant `MCMCState`; `position`/`params`/`potential_energy`/
`gradient`/`inverse_mass_matrix` are the host-side, plot-friendly view of it
(`params` is the *current* draw transformed into constrained model space) — see
[`_current_state`](@ref).

If `checkpoint` is a path, the post-warmup `MCMCState` is written there via
`ProbProg.save_state`, so `sample(...; restart=true)` can skip warmup and continue
sampling. Mid-warmup resume is not supported (warmup is one fused call, and
`MCMCState` does not carry the dual-averaging/Welford accumulators needed to resume
adaptation correctly).
"""
function warmup_single(
        rng, ldf, x0, tpost, sampler::ReactantNUTS;
        callback = default_warmup_callback, checkpoint = nothing
    )

    na = sampler.n_adapts
    na > 0 || throw(ArgumentError("n_adapts must be positive"))

    T = eltype(x0)
    step0 = ConcreteRNumber(T(sampler.init_step_size))
    mass0 = Reactant.to_rarray(ones(T, length(x0)))

    # NOTE: `num_samples = 2` (not 1) sidesteps a Reactant 0.2 mcmc_logpdf MLIR
    # bug where the diagnostics result is declared scalar `tensor<i1>` but the
    # body produces `tensor<1xi1>` — function verification fails. We discard
    # these post-warmup samples anyway, so the choice of >=2 here is a workaround.
    run = function (x0, step0, mass0)
        return ProbProg.mcmc_logpdf(
            rng, ldf, x0, tpost;
            algorithm = :NUTS, num_warmup = na, num_samples = 2,
            step_size = step0, inverse_mass_matrix = mass0,
            max_tree_depth = sampler.max_tree_depth,
            max_delta_energy = sampler.max_delta_energy,
            adapt_step_size = true, adapt_mass_matrix = true,
            strong_zero = sampler.strong_zero
        )
    end
    cfn = Reactant.Compiler.compile(run, (x0, step0, mass0); optimize = :probprog)
    samples, diagnostics, _, state = cfn(x0, step0, mass0)

    cur = _current_state(state, tpost)
    info = (;
        num_warmup = na,
        step_size = cur.step_size,
        position = cur.position, params = cur.params,
        potential_energy = cur.potential_energy, gradient = cur.gradient,
        inverse_mass_matrix = cur.inverse_mass_matrix,
        state, samples = Array(samples), diagnostics = Array(diagnostics),
    )
    history = callback(info)

    isnothing(checkpoint) || ProbProg.save_state(checkpoint, state)
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
            samples, diagnostics, _, state = cfn(state, ns)
        end

        raw = Array(samples)
        chain = [transform(tpost, r) for r in eachrow(raw)]
        numerical_error = .!Array(diagnostics)   # diagnostics: true == NO divergence

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
        HypercubeTransform.inverse(tpost, initial_params)
    end
    return Reactant.to_rarray(x)
end

_default_ldf(x, tpost) = logdensityof(tpost, x)

"""
    sample(rng, post, sampler::ReactantNUTS, nsamples;
           saveto=MemoryStore(), initial_params=nothing, restart=false,
           chunk_size=100, ldf=_default_ldf, host_rng=Random.default_rng(),
           warmup_callback=default_warmup_callback)

Warm up (single fused Stan-windowed adaptation, run internally by the ProbProg NUTS
pass — see [`warmup_single`](@ref)) then draw `nsamples` post-warmup samples from the
Reactant posterior `post`. Structured like the AdvancedHMC extension's `sample`, and
algorithmically identical to AdvancedHMC's NUTS adaptation.

Returns a `NamedTuple` `(; out, state)` where `state` is the final ProbProg
`MCMCState` (held in memory, ready to inspect, plot via [`_current_state`](@ref), or
thread into a follow-up `sample_chunked`) and `out` is the standard Comrade output:

  - `saveto::MemoryStore` -> `out` is a `PosteriorSamples` (chain transformed to
    constrained space; `samplerstats` carries per-sample `numerical_error`;
    warmup/sample history + final state in the metadata).
  - `saveto::DiskStore` -> writes per-chunk `PosteriorSamples` to `saveto.name` in
    Comrade's on-disk layout; `out` is the `Comrade.DiskOutput` handle. Read the
    chain back with `Comrade.load_samples(out)` or `Comrade.load_samples(saveto.name)`.
    A resumable `MCMCState` is checkpointed to `<name>/state.jls` after warmup
    completes and after each sampling chunk.

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

`warmup_callback` runs once, after warmup completes; its `info` carries the
warmup-specific field `num_warmup` alongside the host-side state view (see
[`default_warmup_callback`](@ref) and [`warmup_single`](@ref)).

`restart=true` (only with a `DiskStore`) continues an interrupted run, matching
AdvancedHMC's `restart`. Because warmup is a single fused call, it cannot be resumed
mid-warmup; `restart` therefore requires a completed-warmup `state.jls` checkpoint
(written once warmup finishes) and continues from the sampling phase, appending to
the chain already on disk.

`nsamples` is the TOTAL target chain length, samples already on disk are counted,
and new chunks are numbered *after* them and appended (with `parameters.jls` grown
to the cumulative total). If the requested `nsamples` is already on disk, nothing is
drawn.
"""
function AbstractMCMC.sample(
        rng::Reactant.ReactantRNG, post, sampler::ReactantNUTS,
        nsamples::Int;
        saveto = MemoryStore(), initial_params = nothing, restart::Bool = false,
        chunk_size::Int = 100,
        ldf = _default_ldf, host_rng = Random.default_rng(),
        warmup_callback = default_warmup_callback
    )

    tpost = asflat(post)

    if restart
        saveto isa DiskStore ||
            throw(ArgumentError("restart=true requires saveto::DiskStore"))
        ckpt = joinpath(saveto.name, "state.jls")
        isfile(ckpt) ||
            throw(ArgumentError("cannot restart: no state checkpoint at $ckpt"))
    end

    # Warmup state comes either from a fresh single fused warmup, or (on restart)
    # from the completed-warmup checkpoint — warmup is never resumed mid-flight.
    if !restart
        x0 = _initial_position(host_rng, tpost, initial_params)
        # Checkpoint the post-warmup state when saving to disk so sampling is resumable.
        warmup_ckpt = saveto isa DiskStore ?
            (mkpath(saveto.name); joinpath(saveto.name, "state.jls")) : nothing
        @info "ReactantNUTS warmup" n_adapts = sampler.n_adapts
        state, warmup_history = warmup_single(
            rng, ldf, x0, tpost, sampler;
            callback = warmup_callback, checkpoint = warmup_ckpt
        )
    else
        # Restart path: warmup already completed; load the checkpointed state (the
        # latest sampling state if any chunks ran, else the post-warmup state) and
        # continue sampling.
        @info "ReactantNUTS restart: loading checkpoint, skipping warmup" dir = saveto.name
        state = ProbProg.load_state(joinpath(saveto.name, "state.jls"))
        warmup_history = nothing
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
