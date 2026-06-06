using AbstractMCMC
using Serialization
using Random
using HypercubeTransform
using Printf
using Reactant: ProbProg

"""
    n_adapts(sampler::ReactantNUTS) -> Int

Total number of warmup/adaptation steps. Equals `sampler.n_adapts` for the derived
schedule, or `init_buffer + sum(windows) + term_buffer` when `windows` is overridden.
"""
function n_adapts(s::ReactantNUTS)
    isnothing(s.windows) && return s.n_adapts
    return s.init_buffer + sum(s.windows; init = 0) + s.term_buffer
end

"""
    _stan_windows(num_warmup, init_buffer, term_buffer, base_window)
        -> (init_buffer, windows::Vector{Int}, term_buffer)

Stan's windowed warmup schedule. The middle region (`num_warmup - init_buffer -
term_buffer`) is split into windows that start at `base_window` and double; the final
window is expanded to absorb the remainder so the three pieces sum exactly to
`num_warmup`. If the buffers + one base window don't fit, fall back to Stan's
15%/10% buffer split with a single window.
"""
function _stan_windows(num_warmup::Int, init_buffer::Int, term_buffer::Int, base_window::Int)
    num_warmup > 0 || throw(ArgumentError("n_adapts must be positive"))
    if init_buffer + base_window + term_buffer > num_warmup
        ib = round(Int, 0.15 * num_warmup)
        tb = round(Int, 0.1 * num_warmup)
        w = num_warmup - ib - tb
        return ib, [w], tb
    end
    windowed = num_warmup - init_buffer - term_buffer
    windows = Int[]
    w = base_window
    used = 0
    while used < windowed
        remaining = windowed - used
        # Stan: if this window plus the next (doubled) would overrun the region,
        # extend the current window to the end instead.
        if 3w > remaining
            push!(windows, remaining)
            used = windowed
        else
            push!(windows, w)
            used += w
            w *= 2
        end
    end
    return init_buffer, windows, term_buffer
end

"""
    warmup_schedule(sampler::ReactantNUTS) -> (init_buffer, windows::Vector{Int}, term_buffer)

The concrete warmup schedule this sampler will run: either Stan's derived windows
(from `n_adapts`) or the explicit `windows` override.
"""
function warmup_schedule(s::ReactantNUTS)
    isnothing(s.windows) || return (s.init_buffer, collect(Int, s.windows), s.term_buffer)
    return _stan_windows(s.n_adapts, s.init_buffer, s.term_buffer, s.base_window)
end

# ===========================================================================
# Sample-retention backends (reuse Comrade's MemoryStore / DiskStore configs)
# Both deal in transformed PosteriorSamples. A "sink" is driven by:
#   _open_sink -> create   _write_sink! -> stash a chunk PS   _close_sink -> finalize
# ===========================================================================

# --- MemoryStore: accumulate transformed chains + stats, build one PS at close ---
mutable struct _MemorySink
    chains::Vector{Any}
    nerr::Vector{Any}
end
_open_sink(::MemoryStore, _tpost, _nsamples, _nscans, _stride; append::Bool = false) =
    _MemorySink(Any[], Any[])
function _write_sink!(sink::_MemorySink, ::MemoryStore, chunk)
    push!(sink.chains, chunk.chain)
    push!(sink.nerr, chunk.numerical_error)
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
function _write_sink!(sink::_DiskSink, ::DiskStore, chunk)
    sink.iter += 1
    sink.nsamples_done += length(chunk.chain)
    ps = PosteriorSamples(chunk.chain, (; numerical_error = chunk.numerical_error))
    serialize(
        sink.outbase * Printf.@sprintf("%05d.jls", sink.iter),
        (samples = Comrade.postsamples(ps), stats = Comrade.samplerstats(ps))
    )
    # resumable MCMC state checkpoint (latest wins)
    ProbProg.save_state(joinpath(sink.outdir, "state.jls"), chunk.state)
    # update parameters.jls every chunk so a crash mid-run leaves it current and
    # consistent with the files on disk (AHMC checkpoints its counter every scan).
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
sample_callback = info -> plot(skymodel(post, info.params.sky))
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

Default `warmup_callback`: log one line and return the per-round step size.
"""
function default_warmup_callback(info)
    @info "warmup round $(info.round)/$(info.nrounds)" phase = info.phase num_warmup = info.num_warmup step_size = info.step_size
    return (; info.round, info.phase, info.num_warmup, info.step_size)
end

"""
    default_sample_callback(info) -> NamedTuple

Default `sample_callback`: log one line and return the per-chunk divergence count.
"""
function default_sample_callback(info)
    ndiv = count(info.numerical_error)
    @info "sampling round $(info.round)/$(info.nrounds)" num_samples = info.num_samples step_size = info.step_size n_divergences = ndiv
    return (; info.round, info.num_samples, info.step_size, n_divergences = ndiv)
end

# ===========================================================================
# Engine: windowed warmup + chunked sampling
# ===========================================================================

"""
    _build_schedule(sampler::ReactantNUTS) -> Vector{Tuple{Symbol,Int,Bool,Bool}}

Build the per-round `(phase, num_warmup, adapt_step_size, adapt_mass_matrix)` schedule.
"""
function _build_schedule(sampler::ReactantNUTS)
    init_buffer, windows, term_buffer = warmup_schedule(sampler)
    schedule = Tuple{Symbol, Int, Bool, Bool}[]
    init_buffer > 0 && push!(schedule, (:init, init_buffer, true, false))
    for W in windows
        W > 0 && push!(schedule, (:window, W, true, true))
    end
    term_buffer > 0 && push!(schedule, (:term, term_buffer, true, false))
    isempty(schedule) && throw(ArgumentError("warmup schedule is empty"))
    return schedule
end

_warmup_status_path(state_ckpt::AbstractString) = joinpath(dirname(state_ckpt), "warmup_status.jls")

"""
    warmup_windowed(rng, ldf, x0, tpost, sampler::ReactantNUTS;
                    callback, checkpoint,
                    resume_state, resume_round, resume_schedule, resume_history) -> (state, history)

Run the sampler's Stan-style windowed warmup. `ldf(x, tpost)` is the log-density.
Returns the final `MCMCState` and the vector of `callback` returns (one per round).
The `info` NamedTuple passed to `callback` has fields:
  round, nrounds, phase (:init|:window|:term), num_warmup, adapt_step_size,
  adapt_mass_matrix, step_size, position, params, potential_energy, gradient,
  inverse_mass_matrix, state, samples, diagnostics.
`state` is the raw Reactant `MCMCState`; `position`/`params`/`potential_energy`/
`gradient`/`inverse_mass_matrix` are the host-side, plot-friendly view of it
(`params` is the *current* draw transformed into constrained model space) — see
[`_current_state`](@ref).

If `checkpoint` is a path, the post-round `MCMCState` is written there via
`ProbProg.save_state` after every window AND a `warmup_status.jls` is written next
to it recording `(schedule, round_done, complete, history)`. A crashed warmup can
then be resumed from the next window via `sample(...; restart=true)` — the same
schedule is reused so metric adaptation remains consistent.

Pass `resume_state` (non-`nothing`) + `resume_round` to continue from the round
*after* `resume_round`, using `resume_schedule` (must be identical to the original).
`x0` is unused in resume mode.
"""
function warmup_windowed(
        rng, ldf, x0, tpost, sampler::ReactantNUTS;
        callback = default_warmup_callback, checkpoint = nothing,
        resume_state = nothing, resume_round::Int = 0,
        resume_schedule = nothing, resume_history = Any[]
    )

    # Round schedule: (phase, num_warmup, adapt_step_size, adapt_mass_matrix),
    # derived from n_adapts (Stan-style) unless the sampler overrides `windows`.
    # On resume we use the schedule that was persisted at the start of the run
    # so adaptation windows stay consistent.
    schedule = isnothing(resume_schedule) ? _build_schedule(sampler) : resume_schedule
    nrounds = length(schedule)
    0 <= resume_round <= nrounds ||
        throw(ArgumentError("resume_round=$resume_round out of bounds [0, $nrounds]"))

    # On a fresh run with a checkpoint path, persist the schedule up front so
    # restart knows the original adaptation plan even if no round has finished.
    if !isnothing(checkpoint) && isnothing(resume_state) && resume_round == 0
        serialize(
            _warmup_status_path(checkpoint),
            (; schedule, round_done = 0, complete = false, history = Any[])
        )
    end

    T = eltype(x0)
    step0 = ConcreteRNumber(T(sampler.init_step_size))
    mass0 = Reactant.to_rarray(ones(T, length(x0)))

    # Round 1 uses the rng-form from x0; later rounds use the state-form.
    # NOTE: `num_samples = 2` (not 1) sidesteps a Reactant 0.2 mcmc_logpdf MLIR
    # bug where the diagnostics result is declared scalar `tensor<i1>` but the
    # body produces `tensor<1xi1>` — function verification fails. We discard
    # the warmup samples anyway, so the choice of >=2 here is purely a workaround.
    run_round = function (state_or_nothing, nw, adapt_ss, adapt_mm)
        if state_or_nothing === nothing
            return ProbProg.mcmc_logpdf(
                rng, ldf, x0, tpost;
                algorithm = :NUTS, num_warmup = nw, num_samples = 2,
                step_size = step0, inverse_mass_matrix = mass0,
                max_tree_depth = sampler.max_tree_depth,
                max_delta_energy = sampler.max_delta_energy,
                adapt_step_size = adapt_ss, adapt_mass_matrix = adapt_mm,
                strong_zero = sampler.strong_zero
            )
        else
            return ProbProg.mcmc_logpdf(
                state_or_nothing, ldf, tpost;
                algorithm = :NUTS, num_warmup = nw, num_samples = 2,
                max_tree_depth = sampler.max_tree_depth,
                max_delta_energy = sampler.max_delta_energy,
                adapt_step_size = adapt_ss, adapt_mass_matrix = adapt_mm,
                strong_zero = sampler.strong_zero
            )
        end
    end

    compiled = Dict{Tuple{Bool, Int, Bool, Bool}, Any}()
    state = resume_state
    history = copy(resume_history)

    for round in (resume_round + 1):nrounds
        (phase, nw, adapt_ss, adapt_mm) = schedule[round]
        is_first = state === nothing
        key = (is_first, nw, adapt_ss, adapt_mm)
        cfn = get!(compiled, key) do
            arg = is_first ? nothing : state
            Reactant.Compiler.compile(
                run_round, (arg, nw, adapt_ss, adapt_mm);
                optimize = :probprog
            )
        end
        samples, diagnostics, _, state = is_first ?
            cfn(nothing, nw, adapt_ss, adapt_mm) :
            cfn(state, nw, adapt_ss, adapt_mm)

        cur = _current_state(state, tpost)
        info = (;
            round, nrounds, phase,
            num_warmup = nw, adapt_step_size = adapt_ss, adapt_mass_matrix = adapt_mm,
            step_size = cur.step_size,
            position = cur.position, params = cur.params,
            potential_energy = cur.potential_energy, gradient = cur.gradient,
            inverse_mass_matrix = cur.inverse_mass_matrix,
            state, samples = Array(samples), diagnostics = Array(diagnostics),
        )
        push!(history, callback(info))

        # Checkpoint the post-round state + warmup progress so a crashed warmup
        # can resume from the *next* round with the same schedule.
        if !isnothing(checkpoint)
            ProbProg.save_state(checkpoint, state)
            serialize(
                _warmup_status_path(checkpoint),
                (; schedule, round_done = round, complete = (round == nrounds), history)
            )
        end
    end

    return state, history
end

"""
    sample_chunked(state, ldf, tpost, sampler::ReactantNUTS; num_samples, saveto, chunk_size, callback)
        -> (state, out, history)

Draw `num_samples` post-warmup samples in chunks, threading `state` forward with
adaptation OFF (frozen metric + step size). Each chunk is transformed to constrained
space and routed to `saveto` (`MemoryStore` -> accumulate; `DiskStore` -> write
Comrade layout). `out` is a `PosteriorSamples` (memory) or `Comrade.DiskOutput`
(disk). `callback` runs after EVERY chunk; the `info` it receives has fields:
  round, nrounds, num_samples, step_size, position, params, potential_energy,
  gradient, inverse_mass_matrix, state, samples (raw flat matrix),
  chain (transformed Vector), numerical_error (Vector{Bool}).
`state` is the raw Reactant `MCMCState`; `position`/`params`/`potential_energy`/
`gradient`/`inverse_mass_matrix` are its host-side, plot-friendly view (`params` is
the chunk's last draw in constrained model space) — see [`_current_state`](@ref).
"""
function sample_chunked(
        state, ldf, tpost, sampler::ReactantNUTS;
        num_samples::Int, saveto = MemoryStore(), chunk_size::Int = 100,
        append::Bool = false,
        callback = default_sample_callback, metadata = Dict{Symbol, Any}()
    )

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
        samples, diagnostics, _, state = cfn(state, ns)

        raw = Array(samples)
        chain = [transform(tpost, r) for r in eachrow(raw)]
        numerical_error = .!Array(diagnostics)   # diagnostics: true == NO divergence

        cur = _current_state(state, tpost)
        chunk_nt = (;
            round, nrounds, num_samples = ns,
            step_size = cur.step_size,
            position = cur.position, params = cur.params,
            potential_energy = cur.potential_energy, gradient = cur.gradient,
            inverse_mass_matrix = cur.inverse_mass_matrix,
            state, samples = raw, chain, numerical_error,
        )
        _write_sink!(sink, saveto, chunk_nt)
        push!(history, callback(chunk_nt))
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
           warmup_callback=default_warmup_callback,
           sample_callback=default_sample_callback)

Warm up (Stan-style windowed) then draw `nsamples` post-warmup samples from the
Reactant posterior `post` using ProbProg NUTS. Structured like the AdvancedHMC
extension's `sample`.

Returns a `NamedTuple` `(; out, state)` where `state` is the final ProbProg
`MCMCState` (held in memory, ready to inspect, plot via [`_current_state`](@ref), or
thread into a follow-up `sample_chunked`) and `out` is the standard Comrade output:

  - `saveto::MemoryStore` -> `out` is a `PosteriorSamples` (chain transformed to
    constrained space; `samplerstats` carries per-sample `numerical_error`;
    warmup/sample history + final state in the metadata).
  - `saveto::DiskStore` -> writes per-chunk `PosteriorSamples` to `saveto.name` in
    Comrade's on-disk layout; `out` is the `Comrade.DiskOutput` handle. Read the
    chain back with `Comrade.load_samples(out)` or `Comrade.load_samples(saveto.name)`.
    A resumable `MCMCState` is checkpointed to `<name>/state.jls` after each sampling
    chunk and after each warmup window.

`restart=true` (only with a `DiskStore`) continues an interrupted run, matching
AdvancedHMC's `restart`. The on-disk `warmup_status.jls` records whether warmup
finished and which round was last completed:

  - if warmup did **not** complete, the remaining windows are run using the
    *original* schedule (so step-size + metric adaptation stay consistent),
    picking up from `state.jls`;
  - if warmup **did** complete, it is skipped and we go straight to sampling.

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
        warmup_callback = default_warmup_callback,
        sample_callback = default_sample_callback
    )

    tpost = asflat(post)

    if restart
        saveto isa DiskStore ||
            throw(ArgumentError("restart=true requires saveto::DiskStore"))
        ckpt = joinpath(saveto.name, "state.jls")
        isfile(ckpt) ||
            throw(ArgumentError("cannot restart: no state checkpoint at $ckpt"))
    end

    # Warmup state comes either from the restart checkpoint or a fresh windowed warmup.
    if !restart
        x0 = _initial_position(host_rng, tpost, initial_params)
        # Checkpoint warmup per-window when saving to disk, so it too is resumable.
        warmup_ckpt = saveto isa DiskStore ?
            (mkpath(saveto.name); joinpath(saveto.name, "state.jls")) : nothing
        @info "ReactantNUTS warmup" n_adapts = n_adapts(sampler) schedule = warmup_schedule(sampler)
        state, warmup_history = warmup_windowed(
            rng, ldf, x0, tpost, sampler;
            callback = warmup_callback, checkpoint = warmup_ckpt
        )
    else
        # Restart path: consult warmup_status.jls to decide whether warmup is done
        # or needs to be resumed in its remaining windows.
        warmup_ckpt = joinpath(saveto.name, "state.jls")
        status_path = _warmup_status_path(warmup_ckpt)
        loaded_state = ProbProg.load_state(warmup_ckpt)
        if !isfile(status_path)
            # Legacy / external checkpoint without a status sidecar — assume done.
            @warn "ReactantNUTS restart: no warmup_status.jls found; assuming warmup complete" dir = saveto.name
            state = loaded_state
            warmup_history = nothing
        else
            status = deserialize(status_path)
            if status.complete
                @info "ReactantNUTS restart: warmup complete, skipping" dir = saveto.name round_done = status.round_done
                state = loaded_state
                warmup_history = status.history
            else
                @info "ReactantNUTS restart: resuming warmup" dir = saveto.name resume_from = status.round_done + 1 nrounds = length(status.schedule)
                # x0 is unused in resume mode (state carries position/metric/step), but
                # length/eltype seed the unused-but-constructed step0/mass0 arrays.
                T0 = eltype(Array(loaded_state.position))
                x0_dummy = Reactant.to_rarray(zeros(T0, HypercubeTransform.dimension(tpost)))
                state, warmup_history = warmup_windowed(
                    rng, ldf, x0_dummy, tpost, sampler;
                    callback = warmup_callback, checkpoint = warmup_ckpt,
                    resume_state = loaded_state, resume_round = status.round_done,
                    resume_schedule = status.schedule, resume_history = status.history
                )
            end
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
        num_samples = remaining, saveto, chunk_size, append = restart,
        callback = sample_callback, metadata
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
