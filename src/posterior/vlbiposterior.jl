struct VLBIPosterior{DV, DI, LV, LI, P, MS <: AbstractSkyModel, MI <: AbstractInstrumentModel, ADMode <: Union{Nothing, EnzymeCore.Mode}} <: AbstractVLBIPosterior
    data::DV
    dataimg::DI
    lklhds::LV
    lklhdsimg::LI
    prior::P
    skymodel::MS
    instrumentmodel::MI
    admode::ADMode
end

(post::VLBIPosterior)(θ) = logdensityof(post, θ)
admode(post::VLBIPosterior) = post.admode

@inline function datatype(::Type{<:VLBIPosterior{DV, DI}}) where {DV, DI}
    return datatype(DV, DI)
end


"""
    VLBIPosterior(skymodel::SkyModel, instumentmodel::InstrumentModel, 
                  dataproducts::EHTObservationTable...; 
                  admode=set_runtime_activity(set_strong_zero(Reverse)))

Creates a VLBILikelihood using the `skymodel` its related metadata `skymeta`
and the `instrumentmodel` and its metadata `instumentmeta`. The `model` is a 
function that converts from parameters `θ` to a Comrade
AbstractModel which can be used to compute [`visibilitymap`](@ref) and a set of
`metadata` that is used by `model` to compute the model.

To enable automatic differentiation, the `admode` keyword argument can be set to any `EnzymeCore.Mode` type
of if no AD is desired then `nothing`. The default is `Enzyme.set_runtime_activity(Enzyme.Reverse)`
which should work for essentially every problem. Note that runtime activity does have a perfomance cost,
and as Enzyme and Comrade matures we expect this to not need runtime activity.

!!! note
    When the sky model is defined over a multi-frequency/multi-time image grid (one carrying `Ti`/`Fr`
    dimensions), the visibility data products are automatically reordered at construction so that each
    `Ti`/`Fr` plane's visibilities form a contiguous block matching the image's layout. This improves
    locality for the multidomain Fourier transform and enables sharding across devices. The reordering is
    a permutation, so the posterior is unchanged, but `measurement(post.data[i])` will be in this
    regrouped order rather than the input order.

# Warning

The `model` itself must be a two argument function where the first argument is the set
of model parameters and the second is a container that holds all the additional
information needed to construct the model. An example of this is when the model
needs some precomputed cache to define the model.

# Example
```julia
dlcamp, dcphase = extract_table(obs, LogClosureAmplitude(), ClosurePhases())
array = arrayconfiguration(dlcamp)

function sky(θ, metadata)
    (; r, a) = θ
    m = stretched(ExtendedRing(a), r, r)
    return m
end

skyprior = (r = VLBIUniform(μas2rad(10.0), μas2rad(30.0)), a = VLBIUniform(1.0, 10.0))
g  = imagepixels(μas2rad(100.0), μas2rad(100.0), 256, 256)
skym = SkyModel(sky, skyprior, g)

G = SingleStokesGain(x->exp(complex(x.lg, x.pg)))
intprior = (lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            pg = ArrayPrior(IIDSitePrior(ScanSeg(), DiagVonMises(0.0, inv(π^2))))
            )
intmodel = InstrumentModel(G, intprior, array)

post = VLBIPosterior(skym, intmodel, dlcamp, dcphase)
```
"""
@noinline function VLBIPosterior(
        skymodel::AbstractSkyModel,
        instrumentmodel::AbstractInstrumentModel,
        dataproducts::EHTObservationTable...;
        imgdata = nothing,
        admode = EnzymeCore.set_runtime_activity(EnzymeCore.set_strong_zero(EnzymeCore.Reverse))
    )

    # Sort the visibility data into contiguous Ti/Fr blocks matching the sky image's non-spatial
    # dimensions (see [`_regroup_dataproducts`](@ref)). This improves locality for the per-block
    # multidomain NUFFT and is the precondition for cleanly sharding those blocks across devices.
    dataproducts = _regroup_dataproducts(skymodel, dataproducts)

    array = arrayconfig(dataproducts[begin])
    int, intprior = set_array(instrumentmodel, array)
    sky, skyprior = set_array(skymodel, array)
    total_prior = combine_prior(skyprior, intprior)

    ls = Tuple(map(makelikelihood, dataproducts))
    if !isnothing(imgdata)
        ils = Tuple(map(makelikelihood, imgdata))
    else
        ils = ()
    end

    return VLBIPosterior{
        typeof(dataproducts), typeof(imgdata), typeof(ls), typeof(ils), typeof(total_prior),
        typeof(sky), typeof(int), typeof(admode),
    }(dataproducts, imgdata, ls, ils, total_prior, sky, int, admode)
end

VLBIPosterior(
    skymodel::AbstractSkyModel, dataproducts::EHTObservationTable...;
    admode = EnzymeCore.set_runtime_activity(EnzymeCore.set_strong_zero(EnzymeCore.Reverse)), kwargs...
) =
    VLBIPosterior(skymodel, IdealInstrumentModel(), dataproducts...; admode, kwargs...)


# The non-spatial image axes (`Ti`/`Fr` today, any additional plane axis in future) over which
# `skymodel`'s image is laid out, in regular Julia (column-major) indexing order — fastest image dim
# first, slowest last. These are exactly the grid's `dims[3:end]`, the same non-spatial axes the
# multidomain NUFFT planner iterates, so the two definitions cannot drift apart. A plain `X`/`Y` image
# has none, so grouping is a no-op there.
_grouping_axes(::AbstractSkyModel) = ()
function _grouping_axes(m::Union{SkyModel, FixedSkyModel})
    m.grid isa AbstractRectiGrid || return ()
    return keys(m.grid)[3:end]
end

_isclosuredata(::EHTObservationTable{<:ClosureProducts}) = true
_isclosuredata(::AbstractObservationTable) = false

"""
    _regroup_dataproducts(skymodel, dataproducts)

Reorder each visibility/coherency `dataproduct` so that `Ti` and `Fr` are grouped into contiguous blocks.
The ordering is determined by the image grid in `skymodel` when it has dimensions of `Ti` and `Fr`. 
When the image doesn't have `Ti` or `Fr` dimensions this returns the original array. In the future,
we may regroup the data so that visibilities are group in `Ti` or `Fr`. However, in this case
we will probably rely on the nature cube structure of VLBI data to handle that when you have spectral
windows or multiple scans. 

!!! note
    We currently do not support regrouping of closure products. This will be addressed in the future.
    Visibilities and Coherencies are re-grouped automatically.

!!! note
    A consequence is that `measurement(post.data[i])` is in the regrouped (not input) order. Everything
    derived from the posterior (residuals, simulated data, ...) is self-consistent in that same order.
"""
function _regroup_dataproducts(skymodel::AbstractSkyModel, dataproducts::Tuple)
    axes = _grouping_axes(skymodel)
    (isempty(axes) || isempty(dataproducts)) && return dataproducts
    any(_isclosuredata, dataproducts) && return dataproducts
    # All non-closure products share a single array configuration in Comrade's one-model-visibility
    # architecture, so the permutation derived from the first is applied to every product. Guard that
    # they really are the same length rather than silently dropping rows / throwing a `BoundsError`.
    n = length(dataproducts[begin])
    all(dp -> length(dp) == n, dataproducts) || throw(
        ArgumentError(
            "cannot regroup visibility/coherency data products of differing length " *
                "$(map(length, dataproducts)); they must share one array configuration"
        )
    )
    # The plane ordering lives entirely in `ComradeBase.regroup(domain, grid)`; it derives its axes
    # from the grid's `dims[3:end]` and ranks by grid position, matching the NUFFT's plane enumeration.
    _, perm = ComradeBase.regroup(domain(dataproducts[begin]), skymodel.grid)
    issorted(perm) && return dataproducts
    return map(dp -> dp[perm], dataproducts)
end

function combine_prior(skyprior, instrumentmodelprior)
    return NamedDist((sky = skyprior, instrument = instrumentmodelprior))
end

function combine_prior(skymodel, ::NamedDist{()})
    return NamedDist((; sky = skymodel))
end

function combine_prior(::NamedDist{()}, instrumentmodel)
    return NamedDist((; instrument = instrumentmodel))
end

function combine_prior(::NamedDist{()}, ::NamedDist{()})
    return NamedDist()
end


function Base.show(io::IO, mime::MIME"text/plain", post::VLBIPosterior)
    printstyled(io, "VLBIPosterior"; bold = true, color = :light_magenta)
    println(io)
    show(io, mime, post.skymodel)
    println()
    show(io, mime, post.instrumentmodel)
    println()
    printstyled(io, "Data Products: ", color = :light_green)
    return println(io, map(x -> split(string(datumtype(x)), "{")[1], post.data)...)
    # println(io, "  Prior: ", post.prior)
end


"""
    simulate_observation([rng::Random.AbstractRNG], post::VLBIPosterior, θ; add_thermal_noise=true)

Create a simulated observation using the posterior and its data `post` using the parameter
values `θ`. In Bayesian terminology this is a draw from the posterior predictive distribution.

If `add_thermal_noise` is true then baseline based thermal noise is added. Otherwise, we just
return the model visibilities.
"""
function simulate_observation(rng::Random.AbstractRNG, post::VLBIPosterior, θ; add_thermal_noise = true)
    v0 = last(forward_model(post, θ))
    Σn = _visnoise(first(post.data))
    data = post.data

    # Closures are rather annoying. So rather than directly simulating them under the likelihood
    # we simulate the underlying visibilities thermal error and do regular error propogate.
    if add_thermal_noise
        if eltype(v0) <: Complex
            dv = ComplexVisLikelihood(v0, Σn)
        elseif eltype(v0) <: SMatrix
            dv = CoherencyLikelihood(v0, Σn)
        end
        vis = rand(rng, dv)
    else
        vis = baseimage(v0)
    end

    ms = map(x -> likelihood(x, vis).μ, post.lklhds)
    return ntuple(length(data)) do i
        di = data[i]
        ci = arrayconfig(di)
        return EHTObservationTable{datumtype(di)}(ms[i], noise(di), ci)
    end
end
simulate_observation(post::VLBIPosterior, θ; add_thermal_noise = true) = simulate_observation(Random.default_rng(), post, θ; add_thermal_noise)

_visnoise(d::EHTObservationTable) = noise(d) .^ 2
function _visnoise(d::EHTObservationTable{<:EHTCoherencyDatum})
    n = noise(d)
    return map(x -> x .^ 2, n)
end
_visnoise(d::EHTObservationTable{<:EHTLogClosureAmplitudeDatum}) = getfield(arrayconfig(d), :noise) .^ 2
_visnoise(d::EHTObservationTable{<:EHTClosurePhaseDatum}) = getfield(arrayconfig(d), :noise) .^ 2


"""
    residuals(post::VLBIPosterior, θ)

Compute the residuals for each data product in `post` using the parameter values `θ`.
The resturn objects are `EHTObservationTables`, where the measurements are the residuals.
"""
function residuals(post::VLBIPosterior, p)
    vis = last(forward_model(post, p))
    res = map(x -> residual_data(vis, x), post.data)
    return res
end

"""
    chi2(post::AbstractVLBIPosterior, p; reduce=false)

Returns a tuple of the chi-squared values for each data product in the posterior `post` given the parameters `p`.
Note that the chi-square is not reduced. If you want to reduce it by dividing by the number of data
points then set `reduce=true`. 

Note that reduce doesn't take into account the number of model parameters
which for non-linear models is difficult to define globally.
"""
function chi2(post::AbstractVLBIPosterior, p; reduce = false)
    res = residuals(post, p)
    return map(res) do r
        nd = ndata(r)
        c2 = _chi2(r)
        if reduce
            return c2 / nd
        else
            return c2
        end
    end
end

function _chi2(res::EHTObservationTable)
    c2 = sum(datatable(res)) do d
        r2 = @. abs2(d.measurement / d.noise)
        # Check if residual is NaN which means that the data is missing
        isnan(r2) && return zero(r2)
        return r2
    end
    return c2
end

function _chi2(res::EHTObservationTable{<:EHTCoherencyDatum})
    c2 = sum(datatable(res)) do d
        r2 = @. abs2(d.measurement / d.noise)
        r11 = isnan(r2[1, 1]) ? zero(r2[1, 1]) : r2[1, 1]
        r12 = isnan(r2[1, 2]) ? zero(r2[1, 2]) : r2[1, 2]
        r21 = isnan(r2[2, 1]) ? zero(r2[2, 1]) : r2[2, 1]
        r22 = isnan(r2[2, 2]) ? zero(r2[2, 2]) : r2[2, 2]
        return typeof(r2)(r11, r21, r12, r22)
    end
    return c2

end
