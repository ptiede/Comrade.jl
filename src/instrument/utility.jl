export sites_tuple

"""
    sites_tuple(sites, default; reference=nothing kwargs...)
    sites_tuple(obs::EHTObservation, default; reference=nothing, kwargs...)

Convienence function that will construct a `NamedTuple` of objects
whose names are the `sites` in the observation `obs` or explicitly in the argument
`sites`. The `NamedTuple` will be filled with `default` if no kwargs are defined
otherwise each kwarg (key, value) pair denotes a sites and value pair.

Optionally the user can specify a `reference` sites that will be dropped from the tuple.
This is useful for selecting a reference sites for gain phases

## Examples
```julia-repl
julia> sites = (:AA, :AP, :LM, :PV)
julia> sites_tuple(sites, ScanSeg())
(AA = ScanSeg(), AP = ScanSeg(), LM = ScanSeg(), PV = ScanSeg())
julia> sites_tuple(sites, ScanSeg(); AA = FixedSeg(1.0))
(AA = FixedSeg(1.0), AP = ScanSeg(), LM = ScanSeg(), PV = ScanSeg())
julia> sites_tuple(sites, ScanSeg(); AA = FixedSeg(1.0), PV = TrackSeg())
(AA = FixedSeg(1.0), AP = ScanSeg(), LM = ScanSeg(), PV = TrackSeg())
julia> sites_tuple(sites, Normal(0.0, 0.1); reference=:AA, LM = Normal(0.0, 1.0))
(AP = Normal(0.0, 0.1), LM = Normal(0.0, 1.0), PV = Normal(0.0, 0.1))
```
"""
function sites_tuple(sites::NTuple{N, Symbol}, default; kwargs...) where {N}
    out = map(x->get(kwargs, x, default), st)
    return NamedTuple{st}(out)
end
sites_tuple(dvis::EHTObservation, default; kwargs...) = sites_tuple(Tuple(sites(dvis)), default; kwargs...)
sites_tuple(st::AbstractVector{Symbol}, default; kwargs...) = sites_tuple(Tuple(st), default; kwargs...)
