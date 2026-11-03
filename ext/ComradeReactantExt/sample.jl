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

# Remove, with a warning, any chain stored in `dir`. Stale scan files must not survive a
# fresh open: until the first chunk rewrites the index, `load_samples(dir)` would mix them
# into the new chain.
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

# Open Comrade's on-disk sample layout in `dir`, for the main chain and the warmup log. On
# append, new chunks are numbered after the indexed ones; otherwise any chain in `dir` is
# removed.
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

# Serialize one chunk of constrained draws and their column-oriented stats as the next scan
# file. The index is written separately so callers control the write order.
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
    # The state is saved before the index, so a crash between them leaves the index one chunk
    # behind and a restart overwrites the orphaned scan file instead of duplicating it.
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

# --- Warmup draw log: one draw per warmup chunk ---
# ProbProg records no per-step warmup trace; requesting samples during warmup would add
# non-adapting steps and change the chain. The log is a `_DiskSink` with `stride = 1`, so
# `Comrade.load_samples(dir)` reads it.
_open_warmup_log(dir::String; append::Bool = false) = _open_disk_layout(dir, 1; append)

# `params` is one constrained draw; `stats` is a NamedTuple of length-1 vectors.
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

# Run the `between_chunks` hook. A host position is moved to the device in the compiled
# shape. If the position changed, the cached gradient and potential energy are dropped so the
# next kernel recomputes them.
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

# `Comrade.MoveSet` on a device position: the log density is that of the device posterior,
# the traceable moves of a call run as one program, and other moves propose on the host.
function Comrade._compile_moves(ms::Comrade.MoveSet, tpost, position::Reactant.AbstractConcreteArray)
    data = map(v -> v isa AbstractArray ? Reactant.to_rarray(v) : v, Base.structdiff(ms.ctx, NamedTuple{(:view,)}))
    ctx = Comrade._bind_context(data, tpost)
    zd = Reactant.to_rarray(collect(Float64, vec(Array(position))))
    ld = Reactant.Compiler.compile((tb, x) -> logdensityof(tb, x), (tpost, zd))
    # Distinct arguments: compiling with one array in both places traces them as one input.
    sp = Comrade.space(ctx.view)
    pd = Reactant.Compiler.compile((tb, x, x′) -> Comrade._prior_delta(sp, tb, x, x′), (tpost, zd, copy(zd)))
    tj = Comrade._traced_moves(ms)
    fused = isempty(tj) ? nothing : _compile_fused(ms, tj, ctx, zd)
    return (; data, ld, pd, fused)
end

function Comrade._move_kernels(ms::Comrade.MoveSet, c, tpost, position::Reactant.AbstractConcreteArray)
    ctx = Comrade._bind_context(c.data, tpost)
    host_propose = (j, xh, step) -> Comrade.propose(ms.moves[j], xh, step, ms.ctx)
    function propose(j, x, step)
        x′, r = host_propose(j, Array(x), step)
        return Reactant.to_rarray(x′), r
    end
    function run_fused(x, ℓ, order, steps, logu)
        x, ℓd, logα, accepted = c.fused(
            ctx, x, Reactant.ConcreteRNumber(ℓ), Reactant.to_rarray(order),
            Reactant.to_rarray(steps), Reactant.to_rarray(logu)
        )
        return x, Float64(ℓd), Array(logα), Array(accepted) .> 0
    end
    return (;
        load = p -> Reactant.to_rarray(collect(Float64, vec(Array(p)))),
        logdensity = x -> Float64(c.ld(tpost, x)),
        prior_delta = (x, x′) -> Float64(c.pd(tpost, x, x′)),
        propose,
        host_propose,
        context = ctx,
        fused = run_fused,
        download = Array,
        upload = (xh, x) -> Reactant.to_rarray(xh),
        store = (z, p) -> Reactant.to_rarray(reshape(convert(Array{eltype(p)}, Array(z)), size(p))),
    )
end

# The move flags and step count are runtime arguments: as constants, the compiler unrolls the
# step loop with one copy of the log density per step.
function _compile_fused(ms, tj, ctx, zd)
    moves = Tuple(ms.moves[tj])
    free = Reactant.to_rarray(Float64.(ms.free[tj]))
    n = sum(ms.rounds[tj])
    nsteps = Reactant.ConcreteRNumber(n)
    order = Reactant.to_rarray(ones(Int, n))
    v = Reactant.to_rarray(zeros(n))
    fused = Reactant.Compiler.compile(
        (c, x, ℓ, order, steps, logu, free, nsteps) ->
        _fused_steps(c, moves, free, nsteps, x, ℓ, order, steps, logu),
        (ctx, zd, Reactant.ConcreteRNumber(0.0), order, v, copy(v), free, nsteps)
    )
    return (c, x, ℓ, order, steps, logu) -> fused(c, x, ℓ, order, steps, logu, free, nsteps)
end

# The traceable moves' steps of one call as a traced loop: step `s` is `steps[s]` of
# `moves[order[s]]`, accepted when `logu[s] < log α` (so a NaN rejects). A move with
# `free[c] > 0` is accepted from the prior alone. Returns the final point, its log density,
# and each step's `log α` and acceptance (1 or 0).
function _fused_steps(ctx, moves, free, nsteps, x, ℓ, order, steps, logu)
    logα = zero(steps)
    accepted = zero(steps)
    @trace track_numbers = false for s in 1:nsteps
        x, ℓ, logα, accepted = _fused_step(ctx, moves, free, s, x, ℓ, order, steps, logu, logα, accepted)
    end
    return x, ℓ, logα, accepted
end

function _fused_step(ctx, moves, free, s, x, ℓ, order, steps, logu, logα, accepted)
    c = Comrade._rget(order, s)
    x′, logdet = _propose_at(moves, c, x, Comrade._rget(steps, s), ctx)
    @trace track_numbers = false if Comrade._rget(free, c) > 0
        ℓ′ = ℓ + Comrade._prior_delta(Comrade.space(ctx.view), ctx.view.tbase, x, x′)
    else
        ℓ′ = logdensityof(ctx.view.tbase, x′)
    end
    a = ℓ′ - ℓ + logdet
    accept = Comrade._rget(logu, s) < a
    x = ifelse.(accept, x′, x)
    ℓ = ifelse(accept, ℓ′, ℓ)
    Comrade.ComradeBase.rsetindex!(logα, a, s)
    Comrade.ComradeBase.rsetindex!(accepted, ifelse(accept, 1.0, 0.0), s)
    return x, ℓ, logα, accepted
end

# The proposal of `moves[c]` and its log-determinant: one traced branch per move, in sequence
# rather than nested, so tracing depth does not grow with the number of moves.
function _propose_at(moves, c, x, u, ctx)
    x′ = x
    r = Reactant.promote_to(Reactant.TracedRNumber{Float64}, 0.0)
    for k in eachindex(moves)
        x′, r = _propose_if(moves[k], c == k, x, u, ctx, x′, r)
    end
    return x′, r
end

function _propose_if(m, hit, x, u, ctx, x′, r)
    @trace track_numbers = false if hit
        y, l = _uniform_propose(m, x, u, ctx)
    else
        y, l = x′, r
    end
    return y, l
end

# Every branch of a traced `if` must return the same types.
function _uniform_propose(m, x, u, ctx)
    x′, r = Comrade.propose(m, x, u, ctx)
    return Reactant.promote_to(x, x′), Reactant.promote_to(Reactant.TracedRNumber{Float64}, r)
end

# ===========================================================================
# Default callbacks (called between rounds; return value is collected into history)
# ===========================================================================

"""
    default_warmup_callback(info) -> NamedTuple

Default `warmup_callback` for a `MemoryStore`: log one line per warmup chunk and return its
progress, the current adapted step size, and the current draw (`params`). The return values
are collected into `warmup_history` in `samplerinfo(out)`, which holds the warmup chain.
"""
function default_warmup_callback(info)
    @info "ReactantNUTS warmup" step = info.step total = info.total step_size = info.step_size
    return (; info.step, info.total, info.step_size, info.params)
end

"""
    default_warmup_callback_noparams(info) -> NamedTuple

Default `warmup_callback` for a `DiskStore`: [`default_warmup_callback`](@ref) without
`params` in the return value, since the draws are written to `<name>/warmup` (see
[`warmup_chunked`](@ref)) and keeping them in `warmup_history` would hold the warmup chain in
memory and duplicate it in `metadata.jls`.
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
        # `_infer` returns (trace, diagnostics, log_densities, traced_result, state).
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
`potential_energy`/`gradient`/`inverse_mass_matrix`, see [`_current_state`](@ref)), the
raw `state`, `nsteps` (the NUTS steps since the previous callback), `nuts_seconds` (their
wall time on the device) and `fresh_kernel` (whether one of them ran a just-compiled kernel).

If `checkpoint` is a path, the `MCMCState` is written there with `ProbProg.save_state`
every `checkpoint_every` steps (it persists the adaptation accumulators), and the step
count is recorded to `progress_checkpoint`. Together these let `sample(...; restart=true)`
resume an interrupted warmup from the last checkpoint. To resume, pass the loaded state
as `resume_state` and the recorded step count as `warmup_done`.

`between_chunks(state, tpost, info, host_rng) -> state` runs at the end of every chunk and
may move `state.position` (e.g. a Metropolis–Hastings step); see the `ReactantNUTS`
`sample` method. Here `info` is `(; phase = :warmup, step, total, report)`, where `report`
is true at the stops that checkpoint (or, without checkpointing, at chunk ends) and at the
last step.

If `draws_dir` is a path, the draw at the end of each chunk is appended there as it is
produced, giving an inspectable warmup chain of `n_adapts ÷ chunk` draws:

```julia
warmup = Comrade.load_samples(joinpath("Results", "warmup"))   # PosteriorSamples
Comrade.samplerstats(warmup).step_size                         # adaptation trace
```

Per-draw stats are `step` (global warmup step), `step_size`, and `potential_energy`. A
smaller `chunk` gives a finer log. Passing `resume_state` appends to an existing log rather
than restarting it.
"""
function warmup_chunked(
        rng, ldf, x0, tpost, sampler::ReactantNUTS;
        chunk::Int, callback = default_warmup_callback,
        checkpoint = nothing, progress_checkpoint = nothing,
        draws_dir = nothing,
        resume_state = nothing, warmup_done::Int = 0,
        between_chunks = nothing, host_rng = Random.default_rng(),
        checkpoint_every::Int = chunk,
    )

    na = sampler.n_adapts
    na > 0 || throw(ArgumentError("n_adapts must be positive"))
    chunk > 0 || throw(ArgumentError("warmup chunk length must be positive"))
    checkpoint_every > 0 ||
        throw(ArgumentError("checkpoint_every must be positive, got $checkpoint_every"))

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

    wlog = isnothing(draws_dir) ? nothing :
        _open_warmup_log(draws_dir; append = !isnothing(resume_state))

    history = Any[]
    # One kernel per chunk length and state shape: a between-chunks move empties the gradient
    # slot, so the shapes alternate.
    kernels = Dict{Any, Any}()
    # A kernel runs to the next chunk end or checkpoint, whichever comes first.
    nextmult(k) = (done ÷ k + 1) * k
    since = (; seconds = 0.0, steps = 0, fresh = false)
    while done < na
        stop = min(na, nextmult(chunk))
        isnothing(checkpoint) || (stop = min(stop, nextmult(checkpoint_every)))
        nsteps = stop - done
        key = (nsteps, isnothing(state.gradient), isnothing(state.adaptation), ndims(state.position))
        fresh = !haskey(kernels, key)
        kernel = get!(kernels, key) do
            _compile_warmup_kernel(state, ldf, tpost, nsteps, na, sampler)
        end
        # The kernel runs asynchronously, so synchronize before stopping the timer.
        nuts_seconds = @elapsed begin
            state = kernel(state, ldf, ConcreteRNumber(Int64(done)))
            Reactant.synchronize(state.position)
        end
        done += nsteps
        since = (;
            seconds = since.seconds + nuts_seconds, steps = since.steps + nsteps,
            fresh = since.fresh || fresh,
        )

        # Checkpoint AFTER the steps so a crash resumes from completed work.
        if !isnothing(checkpoint) && (done % checkpoint_every == 0 || done == na)
            ProbProg.save_state(checkpoint, state)
            isnothing(progress_checkpoint) ||
                serialize(progress_checkpoint, (; warmup_done = done, n_adapts = na))
        end
        (done % chunk == 0 || done == na) || continue

        cur = _current_state(state, tpost)
        # Log the draw before the callback, which may throw.
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
            state, nsteps = since.steps, nuts_seconds = since.seconds,
            fresh_kernel = since.fresh,
        )
        push!(history, callback(info))
        since = (; seconds = 0.0, steps = 0, fresh = false)

        report = done == na || isnothing(checkpoint) || done % checkpoint_every == 0
        state = _run_between_chunks(
            between_chunks, state, tpost, (; phase = :warmup, step = done, total = na, report), host_rng
        )
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
`NamedTuple` documented under the `ReactantNUTS` `sample` method (common fields in
[`Comrade.default_disk_callback`](@ref), with the host-side `MCMCState` view from
[`_current_state`](@ref) under `info.extras`).

`between_chunks(state, tpost, info, host_rng) -> state` runs after the callback of every
chunk but the last, with `info = (; phase = :sampling, step, total)`; see the
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

    # Keyed on gradient presence too: a between-chunks move empties the gradient slot.
    compiled = Dict{Tuple{Int, Bool}, Any}()
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
            # (trace, diagnostics, log_densities, traced_result, state)
            samples, diagnostics, log_densities, _, state = cfn(state, ns)
            Reactant.synchronize(state.position)
        end

        raw = Array(samples)
        chain = [transform(tpost, r) for r in eachrow(raw)]
        # `diagnostics` is `ns × 2`. Column 1 is `true` for every draw and carries no
        # information. Column 2 is the divergence signal: 0 for a healthy draw, up to 1 when
        # the energy error exceeds `max_delta_energy`. `numerical_error` must be a
        # `Vector{Bool}` (see `default_disk_callback`), so any nonzero signal is divergent.
        numerical_error = map(>(0), Array(diagnostics)[:, 2])

        # Persist the chunk first (serialization / accumulation), then build the callback
        # `info` and fire the callback — the same ordering the AdvancedHMC path uses.
        cur = _current_state(state, tpost)
        # The step size and the per-chunk wall time repeat over the chunk. `time` excludes
        # compilation but overestimates the first run of a fresh kernel, which pays one-time
        # device setup.
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
                state, samples = raw, fresh_kernel = fresh,
            ),
        )
        push!(history, callback(info))
        ndone += ns

        if round < nrounds
            state = _run_between_chunks(
                between_chunks, state, tpost,
                (; phase = :sampling, step = ndone, total = num_samples), host_rng,
            )
        end
    end

    meta = merge(
        Dict{Symbol, Any}(:nsamples => num_samples, :sample_history => history, :final_state => state),
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
           chunk_size=100, warmup_chunk=0, checkpoint_every=20, ldf=_default_ldf,
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

Returns the same output as the AdvancedHMC extension's `sample`:

  - `saveto::MemoryStore` -> a `PosteriorSamples` (chain transformed to
    constrained space; `samplerstats` carries per-sample `numerical_error`,
    `log_density` (of the flat, unconstrained posterior), `step_size` and `time`
    (wall seconds per draw, averaged over the chunk)). Its `samplerinfo` holds
    `:warmup_history`, `:sample_history` and `:final_state`, the ProbProg `MCMCState`
    after the last draw (see [`_current_state`](@ref) for a host-side view).
  - `saveto::DiskStore` -> writes per-chunk `PosteriorSamples` to `saveto.name` in
    Comrade's on-disk layout and returns the `Comrade.DiskOutput` handle. Read the
    chain back with `Comrade.load_samples(out)` or `Comrade.load_samples(saveto.name)`.
    The same metadata is in `<name>/metadata.jls`. A resumable `MCMCState` is
    checkpointed to `<name>/state.jls` (read it with `ProbProg.load_state`) every
    `checkpoint_every` warmup steps and after each sampling chunk (warmup progress in
    `<name>/warmup_progress.jls`).

## Inspecting warmup

The warmup chain is saved as it runs, one draw per warmup chunk:

  - `saveto::DiskStore` -> the draws are appended to `<name>/warmup` in the same on-disk
    layout as the main chain and load with `Comrade.load_samples(joinpath(name, "warmup"))`.
    Per-draw `samplerstats` carry `step` (global warmup step), `step_size`, and
    `potential_energy`; `Comrade.samplerstats(warmup).step_size` is the step-size
    adaptation trace.
  - `saveto::MemoryStore` -> the draws are the `params` of the entries of
    `samplerinfo(out)[:warmup_history]`.

ProbProg records no per-step warmup trace. Set `warmup_chunk` for a finer log (it defaults
to the sampling chunk size); the chunk length does not otherwise affect the chain.

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
  - `extras.fresh_kernel`        : whether the batch ran a just-compiled kernel, whose first
                                   run pays one-time device setup

`warmup_callback` runs once per warmup chunk; its `info` carries `step`/`total`/
`num_warmup` alongside the host-side state view (see [`default_warmup_callback`](@ref)
and [`warmup_chunked`](@ref)), plus `nsteps` (the NUTS steps of the chunk),
`nuts_seconds` (their wall time on the device, without compilation or moves) and
`fresh_kernel` (whether the chunk ran a just-compiled kernel, whose first run pays one-time
device setup). The default (`nothing`) is [`default_warmup_callback`](@ref) for a
`MemoryStore` and [`default_warmup_callback_noparams`](@ref) for a `DiskStore`, whose draws
are written to `<name>/warmup`. A custom callback's return values replace the entries of
`warmup_history`, so return `info.params` to keep the in-memory warmup chain.

## Moves between chunks

`between_chunks(state, tpost, info, host_rng) -> state` runs on the host between compiled
NUTS chunks and may change `state.position`, e.g. with a Metropolis–Hastings step along a
direction NUTS explores slowly. Composing a posterior-preserving move with the NUTS chunks
preserves the posterior. During warmup it runs at the end of every warmup chunk (after the
callback); during sampling, after the callback of every chunk but the last. `info` is
`(; phase, step, total)`: `phase` is `:warmup` or `:sampling`, and `step` counts the warmup
steps (or post-warmup draws) done so far out of `total`. The constrained point is
`transform(tpost, vec(Array(state.position)))`. During warmup `info` also has `report`, true
at the stops that checkpoint (at chunk ends without checkpointing) and at the last step.
`host_rng` is the `host_rng` keyword.

`state.position` is the latent point of `tpost`; the hook may replace it with a device or
host array of the same length. When the position changes, the next chunk recomputes the
gradient and potential energy; the step size, metric, and adaptation state carry over. The
chain holds only the NUTS draws, and the checkpointed `state.jls`, from which a restart
resumes, is the state before the move. A hook keeps its own acceptance statistics. The
warmup chunk length (by default the `DiskStore` stride) sets how often it runs.

A fresh run (`restart=false`) into a `DiskStore` directory that already holds a chain or
warmup log removes them with a warning before warmup starts.

`restart=true` (only with a `DiskStore`) continues an interrupted run, matching
AdvancedHMC's `restart`. Warmup is checkpointed every `checkpoint_every` steps, so an
interruption *during* warmup resumes from the last checkpoint (the persisted adaptation
state is threaded back in); an interruption during sampling resumes from the last sampling
chunk, appending to the chain already on disk.

`nsamples` is the TOTAL target chain length, samples already on disk are counted,
and new chunks are numbered *after* them and appended (with `parameters.jls` grown
to the cumulative total). If the requested `nsamples` is already on disk, nothing is
drawn.
"""
function AbstractMCMC.sample(
        rng::Reactant.ReactantRNG, post::Comrade.AbstractVLBIPosterior, sampler::ReactantNUTS,
        nsamples::Integer;
        transport_method = nothing,
        saveto = MemoryStore(), initial_params = nothing, restart::Bool = false,
        chunk_size::Integer = 100, warmup_chunk::Integer = 0, checkpoint_every::Integer = 20,
        ldf = _default_ldf, host_rng = Random.default_rng(),
        warmup_callback = nothing, between_chunks = nothing
    )
    nsamples, chunk_size, warmup_chunk, checkpoint_every = Int.((nsamples, chunk_size, warmup_chunk, checkpoint_every))

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
    # The chunk length does not affect the chain: the windowed schedule is anchored to
    # `total_warmup`/`warmup_offset`.
    warmup_chunk = if warmup_chunk > 0
        warmup_chunk
    else
        saveto isa DiskStore ? saveto.stride : chunk_size
    end
    warmup_draws_dir = saveto isa DiskStore ? joinpath(saveto.name, "warmup") : nothing

    if restart
        saveto isa DiskStore ||
            throw(ArgumentError("restart=true requires saveto::DiskStore"))
        isfile(state_ckpt) ||
            throw(ArgumentError("cannot restart: no state checkpoint at $state_ckpt"))
    elseif saveto isa DiskStore
        # Clear before compiling warmup so the warning appears without delay.
        _clear_stale_chain(saveto.name)
        _clear_stale_chain(warmup_draws_dir)
    end

    if !restart
        x0 = _initial_position(host_rng, tpost, initial_params)
        @info "ReactantNUTS warmup" n_adapts = sampler.n_adapts chunk = warmup_chunk
        state, warmup_history = warmup_chunked(
            rng, ldf, x0, tpost, sampler;
            chunk = warmup_chunk, callback = warmup_callback,
            checkpoint = state_ckpt, progress_checkpoint = progress_ckpt,
            draws_dir = warmup_draws_dir, between_chunks, host_rng, checkpoint_every,
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
                between_chunks, host_rng, checkpoint_every,
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
            return deserialize(joinpath(saveto.name, "parameters.jls")).params
        end
        @info "Appending $remaining samples to existing chain of $done" total = nsamples
    end

    metadata = Dict{Symbol, Any}(:sampler => :ReactantNUTS, :warmup_history => warmup_history)
    _, out, _ = sample_chunked(
        state, ldf, tpost, sampler;
        num_samples = remaining, saveto, chunk_size, append = restart, metadata,
        between_chunks, host_rng,
    )
    return out
end

"""
    sample(post, sampler::ReactantNUTS, nsamples; kwargs...)

Convenience method that builds a default `ReactantRNG` and forwards to the
main `sample`.
"""
function AbstractMCMC.sample(post::Comrade.AbstractVLBIPosterior, sampler::ReactantNUTS, nsamples::Integer; kwargs...)
    rng = Reactant.ReactantRNG()
    return AbstractMCMC.sample(rng, post, sampler, nsamples; kwargs...)
end
