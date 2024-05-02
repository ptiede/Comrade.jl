export site_tuple

"""
    site_tuple(sites, default; reference=nothing kwargs...)
    site_tuple(obs::AbstractObservationTable, default; reference=nothing, kwargs...)
    site_tuple(obs::AbstractArrayConfiguration, default; reference=nothing, kwargs...)

Convienence function that will construct a `NamedTuple` of objects
whose names are the `sites` in the observation `obs` or explicitly in the argument
`sites`. The `NamedTuple` will be filled with `default` if no kwargs are defined
otherwise each kwarg (key, value) pair denotes a sites and value pair.

Optionally the user can specify a `reference` sites that will be dropped from the tuple.
This is useful for selecting a reference sites for gain phases

## Examples
```julia-repl
julia> sites = (:AA, :AP, :LM, :PV)
julia> site_tuple(sites, ScanSeg())
(AA = ScanSeg(), AP = ScanSeg(), LM = ScanSeg(), PV = ScanSeg())
julia> site_tuple(sites, ScanSeg(); AA = FixedSeg(1.0))
(AA = FixedSeg(1.0), AP = ScanSeg(), LM = ScanSeg(), PV = ScanSeg())
julia> site_tuple(sites, ScanSeg(); AA = FixedSeg(1.0), PV = TrackSeg())
(AA = FixedSeg(1.0), AP = ScanSeg(), LM = ScanSeg(), PV = TrackSeg())
julia> site_tuple(sites, Normal(0.0, 0.1); reference=:AA, LM = Normal(0.0, 1.0))
(AP = Normal(0.0, 0.1), LM = Normal(0.0, 1.0), PV = Normal(0.0, 0.1))
```
"""
function site_tuple(sites::NTuple{N, Symbol}, default; kwargs...) where {N}
    out = map(x->get(kwargs, x, default), sites)
    return NamedTuple{sites}(out)
end
site_tuple(dvis::AbstractObservationTable, default; kwargs...) = site_tuple(Tuple(sites(dvis)), default; kwargs...)
site_tuple(dvis::AbstractArrayConfiguration, default; kwargs...) = site_tuple(Tuple(sites(dvis)), default; kwargs...)
site_tuple(st::AbstractVector{Symbol}, default; kwargs...) = site_tuple(Tuple(st), default; kwargs...)
