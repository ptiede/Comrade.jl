using AbstractMCMC, PrettyTables
using Serialization

include("optimization.jl")
include("posteriorsamples.jl")


export MemoryStore, DiskStore, load_samples

"""
    MemoryStore

Stored the HMC samplers in memory or ram.
"""
struct MemoryStore end


"""
    DiskStore(;name::String = "Results", stride::Int = 100)

Type that specifies to save the samples results to disk.

# Fields
$(FIELDS)
"""
Base.@kwdef struct DiskStore
    """
    Path of the directory where the results will be saved. If the path does not exist
    it will be automatically created.
    """
    name::String = "Results"
    """
    The output stride, i.e. every `stride` steps the MCMC output will be dumped to disk.
    """
    stride::Int = 100
end
DiskStore(name::String) = DiskStore(name, 100)

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
        indices::Union{Base.Colon, UnitRange, StepRange}=Base.Colon(); table=:both
        )
    @assert (table == :samples || table == :stats || table==:both) "Please select either `samples` or `stats`"
    d = readdir(joinpath(abspath(out.filename), "samples"), join=true)

    if table == :both
        chain = load_samples(out, indices; table=:samples)
        stats = load_samples(out, indices; table=:stats)
        return PosteriorSamples(Comrade.postsamples(chain), stats,
                ; metadata=Dict(:sampler=>:AHMC))
    end


    function tload(f, p)
        s = deserialize(f)
        return getproperty(s, p)
    end

    # load and return the entire table
    if indices == Base.Colon()
        if table == :samples
            return PosteriorSamples(reduce(vcat, tload.(d, :samples)), nothing; metadata=Dict(:sampler=>:AHMC))
        else
            return Comrade.StructArray(reduce(vcat, tload.(d, :stats)))
        end
    end



    # Now get the index of the first file
    ind0 = first(indices)
    (ind0 < 1) && throw(BoundsError(1:out.nsamples, ind0))
    # Now let's find the file
    find0 = ind0÷out.stride + 1
    offset0 = ind0%out.stride # now get the offset
    if offset0 == 0
        find0 = find0 - 1
        offset0 = out.stride
    end

    ind1 = last(indices)
    (ind1 > out.nsamples) && throw(BoundsError(1:out.nsamples, ind1))
    find1 = ind1÷out.stride + 1
    offset1 = ind1%out.stride # now get the offset
    if offset1 == 0
        find1 = find1 - 1
        offset1 = out.stride
    end


    t = reduce(vcat, tload.(d[find0:find1], table))
    out =  t[offset0:step(indices):(out.stride*(find1-find0) + offset1)]
    if table == :samples
        return PosteriorSamples(out, nothing; metadata=Dict(:sampler=>:AHMC))
    else
        return Comrade.StructArray(out)
    end
end

function load_samples(out::String, indices::Union{Base.Colon, UnitRange, StepRange}=Base.Colon(); table=:both)
    @assert isdir(abspath(out)) "$out is not a directory. This isn't where the HMC samples are stored"
    @assert isfile(joinpath(abspath(out), "parameters.jls")) "parameters.jls "
    p = deserialize(joinpath(abspath(out), "parameters.jls")).params
    if p.filename != abspath(out)
        @warn "filename stored in params does not equal what was passed\n"*
                 "we will load the path passed\n  $(out)."
        p = DiskOutput(out, p.nfiles, p.stride, p.nsamples)
    end
    return load_samples(p, indices; table)
end




# include(joinpath(@__DIR__, "fishermatrix.jl"))
