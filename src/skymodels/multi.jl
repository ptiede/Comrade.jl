export MultiSkyModel

"""
    MultiSkyModel(skymodels::NamedTuple)

Create a sky model that is a collection of multiple components. This is useful when
you, e.g., want to decompose the sky into multiple grids, such as for wide-field imaging
where the is a core component and a very far away component. 

# Arguments
 - `skymodels` : A named tuple of sky models, where each model is a [`SkyModel`](@ref).

# Example
```julia
julia> s1 = SkyModel(...)
julia> s2 = SkyModel(...)
julia> stot = MultiSkyModel((core = m1, far = m2))
julia> skymodel(stot, x)
(core = ..., far = ...)
```
"""
struct MultiSkyModel{N, T} <: AbstractSkyModel
    skymodels::NamedTuple{N, T}
end

function ObservedSkyModel(m::MultiSkyModel, arr::AbstractArrayConfiguration)
    skymodels = map(m.skymodels) do sm
        ObservedSkyModel(sm, arr)
    end
    return MultiSkyModel(skymodels)
end

function set_prior(m::MultiSkyModel, array::AbstractArrayConfiguration)
    prs = map(m.skymodels) do sm
        set_prior(sm, array)
    end
    return prs
end

function idealvisibilities(m::MultiSkyModel{N}, x) where {N}
    sm = m.skymodels
    vis = map(N) do n
        Base.@_inline_meta
        @inline idealvisibilities(getproperty(sm, n), (; sky = getproperty(x.sky, n)))
    end
    return reduce(+, vis)
end

function domain(m::MultiSkyModel)
    return map(m.skymodels) do sm
        domain(sm)
    end
end

function skymodel(m::MultiSkyModel{N}, x) where {N}
    sm = m.skymodels
    skyms = map(N) do n
        Base.@_inline_meta
        skymodel(getproperty(sm, n), getproperty(x, n))
    end
    return NamedTuple{N}(skyms)
end

function skymodel(m::MultiSkyModel)
    return map(m.skymodels) do sm
        skymodel(sm)
    end
end
