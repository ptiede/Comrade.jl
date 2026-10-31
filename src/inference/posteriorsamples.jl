export PosteriorSamples, postsamples, samplerstats, samplerinfo, resample_equal


struct PosteriorSamples{T, N, C <: AbstractArray{T, N}, S, M} <: AbstractArray{T, N}
    chain::C
    stats::S
    metadata::M
    function PosteriorSamples(chain, stats::S, metadata::M) where {S, M}
        c = _convert2structarr(chain)
        N = ndims(c)
        T = eltype(c)
        C = typeof(c)
        return new{T, N, C, S, M}(c, stats, metadata)
    end

    function PosteriorSamples(chain::StructArray{T, N}, stats::StructArray, metadata) where {T, N}
        length(chain) != length(stats) && throw(ArgumentError("chain and stats must have the same length"))
        return new{T, N, typeof(chain), typeof(stats), typeof(metadata)}(chain, stats, metadata)
    end

end


# Don't unwrap zero-field types: StructArrays refuses them ("only eltypes with fields are
# supported"), and an empty group (e.g. the sky's `gauss` when no extra Gaussian is fit)
# is fine as a plain column of empty NamedTuples.
_convert2structarr(x::Union{<:NamedTuple, <:AbstractArray}) = StructArray(x, unwrap = (T -> (T <: Union{NamedTuple, Tuple}) && fieldcount(T) > 0))
_convert2structarr(x::StructArray) = x


function PosteriorSamples(chain, stats::AbstractArray; metadata = Dict(:sampler => :unknown))
    ch = _convert2structarr(chain)
    st = StructArray(stats)
    return PosteriorSamples(ch, st; metadata)
end


function PosteriorSamples(chain, stats::NamedTuple; metadata = Dict(:sampler => :unknown))
    ch = _convert2structarr(chain)
    st = StructArray(stats)
    return PosteriorSamples(ch, st; metadata)
end

function PosteriorSamples(chain, stats; metadata = Dict(:sampler => :unknown))
    ch = _convert2structarr(chain)
    return PosteriorSamples(ch, stats, metadata)
end

"""
PosteriorSamples(chain, stats; metadata=Dict(:sampler=>:unknown))

This is the default sampler output from Comrade MCMC extensions. The object contains the
posterior samples, the sampler statistics, and metadata about the sampler used.

Indexing this array behaves like indexing the samples organized as a StructArray. By
default all NamedTuples and Tuples are unwrapped.

The samples of several chains are a draws × chains `PosteriorSamples`, made with `hcat` of
the chains' samples (`reduce(hcat, chains)` for a vector of them). Then `s[i, k]` is draw `i`
of chain `k`, `s[:, k]` is chain `k`, each parameter (e.g. `s.sky.σ`) is a draws × chains
array, and linear indexing, iteration and [`Comrade.rmap`](@ref) run over the draws of all
chains. With MCMCDiagnosticTools loaded, `ess`, `rhat` and `mcse` of a `PosteriorSamples`
give the diagnostic of every parameter, nested like the output of `rmap`.

To access the samples use [`postsamples`](@ref) and to access the sampler statistics use [`samplerstats`](@ref),
or the acesss sampler specific information use [`samplerinfo`](@ref).

To recursively map a function over the samples the unexported [`Comrade.rmap`](@ref).
"""
function PosteriorSamples(chain::StructArray{T, N}, stats::StructArray; metadata = Dict(:sampler => :unknown)) where {T, N}
    length(chain) != length(stats) && throw(ArgumentError("chain and stats must have the same length"))
    return PosteriorSamples(chain, stats, metadata)
end


function Base.parent(s::PosteriorSamples)
    return postsamples(s)
end

function Base.size(s::PosteriorSamples)
    return size(parent(s))
end

Base.getindex(s::PosteriorSamples, i::Integer) = getindex(postsamples(s), i)
Base.getindex(s::PosteriorSamples, i::Integer, j::Integer, I::Integer...) = getindex(postsamples(s), i, j, I...)

# Base's fallbacks for selections and concatenation would copy the samples but leave the
# sampler statistics uninitialized.
function Base.getindex(s::PosteriorSamples, I::Union{Colon, AbstractArray, Integer}...)
    return PosteriorSamples(postsamples(s)[I...], _select(samplerstats(s), I...), samplerinfo(s))
end

_select(::Nothing, I...) = nothing
_select(stats, I...) = stats[I...]

function Base.setindex!(s::PosteriorSamples, v, i::Int)
    return setindex!(postsamples(s), v, i)
end

const VecOrMatSamples = Union{PosteriorSamples{<:Any, 1}, PosteriorSamples{<:Any, 2}}
Base.hcat(s::VecOrMatSamples, ss::VecOrMatSamples...) = _cat(hcat, s, ss...)
Base.vcat(s::VecOrMatSamples, ss::VecOrMatSamples...) = _cat(vcat, s, ss...)
Base.reduce(::typeof(hcat), ss::AbstractVector{<:VecOrMatSamples}) = hcat(ss...)
Base.reduce(::typeof(vcat), ss::AbstractVector{<:VecOrMatSamples}) = vcat(ss...)

# The result carries the metadata of the first samples.
function _cat(f, ss::PosteriorSamples...)
    stats = map(samplerstats, ss)
    all(isnothing, stats) || !any(isnothing, stats) ||
        throw(ArgumentError("cannot concatenate samples with and without sampler statistics"))
    st = all(isnothing, stats) ? nothing : f(stats...)
    return PosteriorSamples(f(map(postsamples, ss)...), st, samplerinfo(first(ss)))
end


function Base.similar(s::PosteriorSamples, ::Type{S}, dims::Dims) where {S}
    isnothing(samplerstats(s)) && return PosteriorSamples(similar(postsamples(s), S, dims), nothing; metadata = samplerinfo(s))
    return PosteriorSamples(similar(postsamples(s), S, dims), similar(samplerstats(s), dims); metadata = samplerinfo(s))
end

Base.IndexStyle(::Type{<:PosteriorSamples{T, N, C}}) where {T, N, C} = Base.IndexStyle(C)

rmap(f, x) = f(x)

"""
    rmap(f, x::PosteriorSamples)

Recursively map a function `f` over the elements of `x`. For instance to compute the mean
of all fields you can do `Comrade.rmap(mean, chain)`
"""
rmap(f, x::PosteriorSamples) = rmap(f, postsamples(x))

function rmap(f, t::Tuple)
    return map(x -> rmap(f, x), t)
end

function rmap(f, nt::NamedTuple{N, T}) where {N, T}
    return NamedTuple{N}(map(x -> rmap(f, x), values(nt)))
end

function rmap(f, x::StructArray)
    return NamedTuple{propertynames(x)}(map(x -> rmap(f, x), StructArrays.components(x)))
end

function rmap(f, x::StructArray{<:Tuple})
    return map(x -> rmap(f, x), StructArrays.components(x))
end

# A parameter group with no parameters (an optional model component that is switched off, e.g.
# a sky model's `gauss` when no extra Gaussian is fitted) has nothing to reduce. StructArrays
# cannot split it into component arrays, so it arrives here as a plain `Vector{@NamedTuple{}}`
# and the generic `rmap(f, x) = f(x)` would call e.g. `mean` on it — which fails in `sum(...)/n`
# with `no method matching /(::@NamedTuple{}, ::Int64)`. Skip `f` and return the empty group.
# (The `StructArray` methods resolve what would otherwise be ambiguities with the ones above.)
# Both empty-group shapes are covered — empty NamedTuples AND empty Tuples — mirroring the
# `fieldcount(T) > 0` unwrap predicate in `_convert2structarr`, which deliberately admits both.
rmap(::Any, ::AbstractArray{<:NamedTuple{(), Tuple{}}}) = NamedTuple()
rmap(::Any, ::StructArray{<:NamedTuple{(), Tuple{}}}) = NamedTuple()
rmap(::Any, ::AbstractArray{Tuple{}}) = ()
rmap(::Any, ::StructArray{Tuple{}}) = ()


function Base.show(io::IO, ::MIME"text/plain", s::PosteriorSamples)
    println(io, "PosteriorSamples")
    println(io, "  Samples size: ", ndims(s) == 2 ? "$(size(s, 1)) draws × $(size(s, 2)) chains" : size(s))
    println(io, "  sampler used: ", get(samplerinfo(s), :sampler, "unknown"))
    ct = postsamples(s)
    pretty_table(
        io, [rmap(mean, ct)]; title = "Mean", alignment = :l
    )
    return pretty_table(
        io, [rmap(std, ct)]; title = "Std. Dev.", alignment = :l
    )

end


function Base.propertynames(s::PosteriorSamples)
    return propertynames(postsamples(s))
end

function Base.getproperty(s::PosteriorSamples, p::Symbol)
    return getproperty(postsamples(s), p)
end

"""
    postsamples(s::PosteriorSamples)

Return the samples from the PosteriorSamples object `s`
"""
function postsamples(s::PosteriorSamples)
    return getfield(s, :chain)
end

"""
    samplerstats(s::PosteriorSamples)

Return the sampler statistics from the PosteriorSamples object `s`
"""
function samplerstats(s::PosteriorSamples)
    return getfield(s, :stats)
end


"""
    samplerinfo(s::PosteriorSamples)

Return the metadata from the PosteriorSamples object `s`.
"""
function samplerinfo(s::PosteriorSamples)
    return getfield(s, :metadata)
end


"""
    resample_equal(post::PosteriorSamples, nsamples::Int)

Resample the posterior samples so you have `nsamples` of equal weight. In order for this method to
be applicable a `:weights` field must be present in the sampler statistics and the weight
must correspond to the probability weights of the samples.
"""
function resample_equal(post::PosteriorSamples, n::Int)
    !(:weights ∈ propertynames(samplerstats(post))) && throw(ArgumentError("Weights not in chain stats, cannot resample"))
    echain = sample(post, AbstractMCMC.StatsBase.Weights(samplerstats(post).weights), n)
    metadata = samplerinfo(post)
    metadata[:sampler] = :NestedSamplers_resampled
    return PosteriorSamples(echain, nothing; metadata)
end
