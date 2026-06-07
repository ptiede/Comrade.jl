using AbstractMCMC, PrettyTables
using Serialization

include("optimization.jl")
include("posteriorsamples.jl")
include("reactant.jl")


export MemoryStore, DiskStore, default_disk_callback, load_samples, ReactantNUTS

"""
    MemoryStore

Stored the HMC samplers in memory or ram.
"""
struct MemoryStore end


"""
    default_disk_callback(info) -> NamedTuple

Default per-batch callback used by [`DiskStore`](@ref). It logs one summary line for the
batch just written to disk and returns a small `NamedTuple` (collected into the sampler's
`sample_history` metadata).

`info` is a `NamedTuple` describing the batch. The fields **guaranteed by every sampler
backend** (AdvancedHMC's `sample` and `ReactantNUTS`) — i.e. the ones this default callback
uses — are:

  - `round::Int`                    : batch index (1-based)
  - `nrounds::Int`                  : total number of batches
  - `num_samples::Int`              : draws in this batch
  - `time::Real`                    : wall-clock seconds spent drawing this batch
  - `numerical_error::Vector{Bool}` : per-draw divergence flags
  - `params`                        : last draw, as the constrained `(; sky[, instrument])` `NamedTuple`
  - `step_size::Real`               : current leapfrog step size
  - `extras::NamedTuple`            : backend-specific data — present on every backend, but
                                      its CONTENTS are NOT guaranteed to match across backends

`info.extras` is where each backend exposes whatever it can cheaply provide (AdvancedHMC's
`samplerstats`, ReactantNUTS's host-side `MCMCState` view, ...). A callback that reaches into
`extras` is therefore backend-specific; the fields are documented on the corresponding
`sample` method (AdvancedHMC's `sample` and `ReactantNUTS`'s `sample`).

A custom callback can do anything with `info`, e.g. plot the current draw:

```julia
DiskStore(; name = "Results", callback = info -> display(imageviz(intensitymap(skymodel(post, info.params.sky), g))))
```
"""
function default_disk_callback(info)
    ndiv = count(info.numerical_error)
    @info "sampling batch $(info.round)/$(info.nrounds)" num_samples = info.num_samples time = info.time step_size = info.step_size n_divergences = "$(ndiv)/$(info.num_samples)"
    return (; info.round, info.num_samples, info.step_size, n_divergences = ndiv)
end


"""
    DiskStore(;name::String = "Results", stride::Int = 100, callback = default_disk_callback)

Type that specifies to save the samples results to disk.

# Fields
$(FIELDS)
"""
Base.@kwdef struct DiskStore{C}
    """
    Path of the directory where the results will be saved. If the path does not exist
    it will be automatically created.
    """
    name::String = "Results"
    """
    The output stride, i.e. every `stride` steps the MCMC output will be dumped to disk.
    """
    stride::Int = 100
    """
    A callback `info -> ...` run after every batch is written to disk. Shared by the
    AdvancedHMC and `ReactantNUTS` disk-sampling paths; see [`default_disk_callback`](@ref)
    for the fields of `info`. The default logs a one-line summary per batch.
    """
    callback::C = default_disk_callback
end
# Positional convenience constructors that preserve the pre-`callback` API. Adding the
# `callback` field made the `@kwdef`-generated positional constructor 3-arg, so restore the
# 1- and 2-arg forms explicitly (they default the callback).
DiskStore(name::String) = DiskStore(name, 100, default_disk_callback)
DiskStore(name::String, stride::Int) = DiskStore(name, stride, default_disk_callback)

struct DiskOutput
    filename::String
    nfiles::Int
    stride::Int
    nsamples::Int
end

"""
    load_samples(out::DiskOutput, indices::Union{Base.Colon, UnitRange, StepRange}=Base.Colon(); table=:samples)
    load_samples(out::String, indices::Union{Base.Colon, UnitRange, StepRange}=Base.Colon(); table=:samples)

The the results from a HMC run saved to disk. To read in the output the user can either
pass the resulting `out` object, or the path to the directory that the results were saved,
i.e. the path specified in [`DiskStore`](@ref).

# Arguments
  - `out::Union{String, DiskOutput}`: If `out` is a string is must point to the direct that the `DiskStore`
     pointed to. Otherwise it is what is directly returned from sample.
  - `indices`: The indices of the that you want to load into memory. The default is to load the entire table.


# Keyword Arguments
  - `table`: A string specifying the table you wish to read in. There are two options: :samples which
     corresponds to the actual MCMC chain, and `stats` which corresponds to additional information
     about the sampler, e.g., the log density of each sample and tree statistics.
"""
function load_samples(
        out::DiskOutput,
        indices::Union{Base.Colon, UnitRange, StepRange} = Base.Colon(); table = :both
    )
    @assert (table == :samples || table == :stats || table == :both) "Please select either `samples` or `stats`"
    d = readdir(joinpath(abspath(out.filename), "samples"), join = true)

    if table == :both
        chain = load_samples(out, indices; table = :samples)
        stats = load_samples(out, indices; table = :stats)
        return PosteriorSamples(
            Comrade.postsamples(chain), stats,
            ; metadata = Dict(:sampler => :AHMC)
        )
    end


    function tload(f, p)
        s = deserialize(f)
        return getproperty(s, p)
    end

    # load and return the entire table
    if indices == Base.Colon()
        if table == :samples
            return PosteriorSamples(reduce(vcat, tload.(d, :samples)), nothing; metadata = Dict(:sampler => :AHMC))
        else
            return Comrade.StructArray(reduce(vcat, tload.(d, :stats)))
        end
    end


    # Now get the index of the first file
    ind0 = first(indices)
    (ind0 < 1) && throw(BoundsError(1:out.nsamples, ind0))
    # Now let's find the file
    find0 = ind0 ÷ out.stride + 1
    offset0 = ind0 % out.stride # now get the offset
    if offset0 == 0
        find0 = find0 - 1
        offset0 = out.stride
    end

    ind1 = last(indices)
    (ind1 > out.nsamples) && throw(BoundsError(1:out.nsamples, ind1))
    find1 = ind1 ÷ out.stride + 1
    offset1 = ind1 % out.stride # now get the offset
    if offset1 == 0
        find1 = find1 - 1
        offset1 = out.stride
    end


    t = reduce(vcat, tload.(d[find0:find1], table))
    out = t[offset0:step(indices):(out.stride * (find1 - find0) + offset1)]
    if table == :samples
        return PosteriorSamples(out, nothing; metadata = Dict(:sampler => :AHMC))
    else
        return Comrade.StructArray(out)
    end
end

function load_samples(out::String, indices::Union{Base.Colon, UnitRange, StepRange} = Base.Colon(); table = :both)
    @assert isdir(abspath(out)) "$out is not a directory. This isn't where the HMC samples are stored"
    @assert isfile(joinpath(abspath(out), "parameters.jls")) "parameters.jls "
    p = deserialize(joinpath(abspath(out), "parameters.jls")).params
    if p.filename != abspath(out)
        @warn "filename stored in params does not equal what was passed\n" *
            "we will load the path passed\n  $(out)."
        p = DiskOutput(out, p.nfiles, p.stride, p.nsamples)
    end
    return load_samples(p, indices; table)
end


# include(joinpath(@__DIR__, "fishermatrix.jl"))
