export PosteriorSamples, postsamples, samplerstats, samplerinfo, resample_equal


struct PosteriorSamples{T, N, C<:AbstractArray{T,N}, S, M} <: AbstractArray{T,N}
    chain::C
    stats::S
    metadata::M
    function PosteriorSamples(chain, stats::S, metadata::M) where {S,M}
        c = _convert2structarr(chain)
        N = ndims(c)
        T = eltype(c)
        C = typeof(c)
        new{T,N,C,S, M}(c, stats, metadata)
    end
    function PosteriorSamples(chain::StructArray{T,N}, stats::StructArray, metadata) where {T,N}
        length(chain) != length(stats) && throw(ArgumentError("chain and stats must have the same length"))
        new{T, N, typeof(chain), typeof(stats), typeof(metadata)}(chain, stats, metadata)
    end

end


_convert2structarr(x::Union{<:NamedTuple, <:AbstractArray}) = StructArray(x, unwrap=(T->(T<:Union{NamedTuple, Tuple})))
_convert2structarr(x::StructArray) = x


function PosteriorSamples(chain, stats::AbstractArray; metadata=Dict(:sampler=>:unknown))
    ch = _convert2structarr(chain)
    st = StructArray(stats)
    return PosteriorSamples(ch, st; metadata)
end


function PosteriorSamples(chain, stats::NamedTuple; metadata=Dict(:sampler=>:unknown))
    ch = _convert2structarr(chain)
    st = StructArray(stats)
    return PosteriorSamples(ch, st; metadata)
end

function PosteriorSamples(chain, stats; metadata=Dict(:sampler=>:unknown))
    ch = _convert2structarr(chain)
    return PosteriorSamples(ch, stats, metadata)
end

function PosteriorSamples(chain::StructArray{T,N}, stats::StructArray; metadata=Dict(:sampler=>:unknown)) where {T,N}
    length(chain) != length(stats) && throw(ArgumentError("chain and stats must have the same length"))
    return PosteriorSamples(chain, stats, metadata)
end




function Base.parent(s::PosteriorSamples)
    return postsamples(s)
end

function Base.size(s::PosteriorSamples)
    return size(parent(s))
end

function Base.getindex(s::PosteriorSamples, i::Int)
    return getindex(postsamples(s), i)
end

# function Base.getindex(s::PosteriorSamples{T,N}, i::Vararg{Int, N}) where {T,N}
#     return getindex(postsamples(s), i)
# end


function Base.setindex!(s::PosteriorSamples, v, i::Int)
    return setindex!(postsamples(s), v, i)
end

# function Base.setindex!(s::PosteriorSamples{T,N}, v, i::Vararg{Int,N}) where {T, N}
#     return setindex!(postsamples(s), v, i)
# end


function Base.similar(s::PosteriorSamples, ::Type{S}, dims::Dims) where {S}
    isnothing(samplerstats(s)) && return PosteriorSamples(similar(postsamples(s), S, dims), nothing; metadata= samplerinfo(s))
    return PosteriorSamples(similar(postsamples(s), S, dims), similar(samplerstats(s), dims); metadata=samplerinfo(s))
end

Base.IndexStyle(::Type{<:PosteriorSamples{T,N,C}}) where {T,N,C} = Base.IndexStyle(C)

rmap(f, x) = f(x)
rmap(f, x::PosteriorSamples) = rmap(f, postsamples(x))

function rmap(f, t::Tuple)
    map(x -> rmap(f,x), t)
end

function rmap(f, nt::NamedTuple{N,T}) where {N,T}
    NamedTuple{N}(map(x -> rmap(f,x), values(nt)))
end

function rmap(f, x::StructArray)
    return NamedTuple{propertynames(x)}(map(x->rmap(f, x), StructArrays.components(x)))
end

function rmap(f, x::StructArray{<:Tuple})
    return map(x->rmap(f, x), StructArrays.components(x))
end


function Base.show(io::IO, ::MIME"text/plain", s::PosteriorSamples)
    println(io, "PosteriorSamples")
    println(io, "  Samples size: $(size(s))")
    println(io, "  sampler used: ", get(samplerinfo(s), :sampler, "unknown"))
    ct = postsamples(s)
    pretty_table(io, [rmap(mean, ct)]; title="Mean", alignment=:l
                )
    pretty_table(io, [rmap(std, ct)]; title="Std. Dev.", alignment=:l
                )

end


function Base.propertynames(s::PosteriorSamples)
    return propertynames(postsamples(s))
end

function Base.getproperty(s::PosteriorSamples, p::Symbol)
    return getproperty(postsamples(s), p)
end

function postsamples(s::PosteriorSamples)
    return getfield(s, :chain)
end

function samplerstats(s::PosteriorSamples)
    return getfield(s, :stats)
end

function samplerinfo(s::PosteriorSamples)
    return getfield(s, :metadata)
end

function resample_equal(post::PosteriorSamples, n::Int)
    !(:weights ∈ propertynames(samplerstats(post))) && throw(ArgumentError("Weights not in chain stats, cannot resample"))
    echain = sample(post, AbstractMCMC.StatsBase.Weights(samplerstats(post).weights), n)
    metadata = samplerinfo(post)
    metadata[:sampler] = :NestedSamplers_resampled
    return PosteriorSamples(echain, nothing; metadata)
end
