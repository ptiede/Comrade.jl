function Comrade.prepare_device(post::VLBIPosterior, ex::ReactantEx)
    for p in propertynames(post)
        p == :prior && continue
        post = Accessors.set(post, PropertyLens(p), Comrade.prepare_device(getproperty(post, p), ex))
    end
    return post
end

function Comrade.prepare_device(m, ex::ReactantEx)
    return Reactant.to_rarray(m)
end

# The posterior's data products and likelihoods are stored as tuples; move each onto the device
# individually so per-element rules (e.g. the sharded `ConditionedLikelihood`/`EHTObservationTable`
# below) get dispatched.
Comrade.prepare_device(t::Tuple, ex::ReactantEx) = map(Base.Fix2(Comrade.prepare_device, ex), t)

# Keep the stored data products coherent with their (sharded) likelihoods. A visibility / coherency
# table's likelihood shards `obs`/`kernel.S`, so shard the *same* measurement and noise vectors here
# onto the same mesh blocks; otherwise `post.data[i]` and `post.lklhds[i]` would hold the measurement
# in two different device layouts and code that re-derives from `post.data` (residuals, simulated data)
# would mix sharded and replicated arrays. Every other product (closures, amplitudes, image-domain)
# is moved unsharded via the generic fallback, matching its unsharded likelihood.
Comrade.prepare_device(d::Comrade.EHTObservationTable, ex::ReactantEx) = _prepare_dataproduct(d, ex)

_prepare_dataproduct(d::Comrade.EHTObservationTable, ::ReactantEx) = Reactant.to_rarray(d)

function _prepare_dataproduct(
        d::Comrade.EHTObservationTable{<:Union{Comrade.EHTVisibilityDatum, Comrade.EHTCoherencyDatum}},
        ex::ReactantEx
    )
    meas = ComradeBase.to_sharded(Comrade.measurement(d), ex)
    noise = ComradeBase.to_sharded(Comrade.noise(d), ex)
    # The array configuration's coordinates are host-iterated by the NUFFT planner (like the grid), so
    # leave them replicated exactly as the previous whole-table `to_rarray` did — only the measurement
    # and noise, which the likelihood shards, need to match the sharded model layout.
    config = Reactant.to_rarray(Comrade.arrayconfig(d))
    return Comrade.EHTObservationTable{Comrade.datumtype(d)}(meas, noise, config)
end

function Comrade.prepare_device(m::Comrade.ObservedSkyModel, ex::ReactantEx)
    grid = Comrade.prepare_device(m.grid, ex)
    return Comrade.ObservedSkyModel(Comrade.prepare_device(m.f, ex), grid, Reactant.to_rarray(m.metadata))
end

function Comrade.prepare_device(m::Comrade.ObservedInstrumentModel, ex::ReactantEx)
    return Comrade.ObservedInstrumentModel(Reactant.to_rarray(m.instrument), m.refbasis, m.bsitelookup)
end

function Comrade.prepare_device(grid::VLBISkyModels.FourierDualDomain, ex::ReactantEx)
    gimr = @jit identity(grid.imgdomain)
    # Shard the (already `regroup`ed, Step-1) visibility coordinate arrays across the mesh declared by
    # `ex`. With no `ShardSpec` this is just `Reactant.to_rarray`, so the unsharded path is unchanged.
    # The uniform split is only a *layout hint*: XLA's sharding propagation inserts whatever resharding
    # the per-plane NUFFT needs, so results stay correct even when the regrouped Ti/Fr blocks are ragged
    # (the common case for real data) and do not line up with the mesh. Block alignment is therefore a
    # best-effort locality optimization, not a correctness requirement — measure before hand-tuning it.
    guvr = ComradeBase.to_sharded(grid.visdomain, ex)
    algr = _device_algorithm(grid.algorithm, grid.imgdomain)
    return VLBISkyModels.FourierDualDomain(gimr, guvr, algr)
end

# Pick the on-device Fourier algorithm. A CPU NUFFT plan can't run under Reactant, so a numeric
# algorithm is swapped for its `ReactantNUFFTAlg` equivalent. `AnalyticAlg` (and any algorithm that
# builds no NUFFT plan) is left untouched: forcing it to `ReactantNUFFTAlg` would call `plan_indices`,
# which iterates host-side over the (now device-resident) visibility coordinates and errors — and it
# defeats the whole point of `AnalyticAlg`, which is to skip FT planning for analytic models.
_device_algorithm(alg::Union{VLBISkyModels.ReactantNUFFTAlg, Comrade.AnalyticAlg}, imgdomain) = alg
_device_algorithm(alg, imgdomain) = VLBISkyModels.ReactantNUFFTAlg(eltype(imgdomain))


# Visibility / coherency likelihoods store the data (`obs`) and noise (`kernel.S`) as flat
# per-visibility vectors aligned 1:1 with the regrouped visibility domain, so we shard them onto the
# same mesh blocks as the model visibilities (`kernel.L` is a scalar log-norm, left replicated). Every
# other likelihood — closures, amplitudes, image-domain — couples its data through a design matrix /
# reduction, so it is moved to the device unsharded via the generic fallback.
function Comrade.prepare_device(d::Comrade.ConditionedLikelihood, ex::ReactantEx)
    return _prepare_likelihood(d.kernel, d, ex)
end

_prepare_likelihood(@nospecialize(kernel), d::Comrade.ConditionedLikelihood, ex::ReactantEx) =
    Reactant.to_rarray(d)

function _prepare_likelihood(
        kernel::Union{Comrade._Visibility, Comrade._Coherency},
        d::Comrade.ConditionedLikelihood, ex::ReactantEx
    )
    knew = Accessors.@set kernel.S = ComradeBase.to_sharded(kernel.S, ex)
    obs = ComradeBase.to_sharded(d.obs, ex)
    return Comrade.ConditionedLikelihood(knew, obs)
end
