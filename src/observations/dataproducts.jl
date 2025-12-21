export extract_table,
    ClosurePhases, LogClosureAmplitudes,
    VisibilityAmplitudes, Visibilities, Coherencies


# Traits that decide the domain of the data. We include this to prevent additional computation if
# e.g., on visibility data is considered.
struct NoData end
struct DualData end
struct VisData end
struct ImgData end

datatype(::Type, ::Nothing) = VisData()
datatype(::Nothing, ::Type) = ImgData()
datatype(::Nothing, ::Nothing) = throw(ArgumentError("No data in uv plane or image plane is provided"))
datatype(::Type, ::Type) = DualData()


abstract type VLBIDataProducts{K} end

keywords(d::VLBIDataProducts) = d.keywords

"""
    ClosuresPhases(;kwargs...)

Type to specify to extract the closure phase table in the [`extract_table`](@ref) function.
Optional keywords are passed through `extract_table` to specify additional option.

# Special keywords for eht-imaging with Pyehtim.jl
For a list of potential keyword arguments see `eht-imaging` and `add_cphase` command for obsdata.
In addition note we have changed the following:
 - count: How the closures are formed, the available options are "min-correct", "min", "max"

## Warning

The `count` keyword argument is treated specially in `Comrade`. The default option
is "min-correct" and should almost always be used.
This option construct a minimal set of closure phases that is valid even when
the array isn't fully connected. For testing and legacy reasons we `ehtim` other count
options are also included. However, the current `ehtim` count="min" option is broken
and does construct proper minimal sets of closure quantities if the array isn't fully connected.

"""
struct ClosurePhases{K} <: VLBIDataProducts{K}
    keywords::K
end

"""
    LogClosureAmplitudes(;kwargs...)

Type to specify to extract the log closure amplitudes table in the [`extract_table`](@ref) function.
Optional keywords are passed through `extract_table` to specify additional option.

# Special keywords for eht-imaging with Pyehtim.jl
For a list of potential keyword arguments see `eht-imaging` and `add_cphase` command for obsdata.
In addition note we have changed the following:
 - count: How the closures are formed, the available options are "min-correct", "min", "max"

Returns an EHTObservation with log-closure amp. datums

## Warning
The `count` keyword argument is treated specially in `Comrade`. The default option
is "min-correct" and should almost always be used.
This option construct a minimal set of closure phases that is valid even when
the array isn't fully connected. For testing and legacy reasons we `ehtim` other count
options are also included. However, the current `ehtim` count="min" option is broken
and does construct proper minimal sets of closure quantities if the array isn't fully connected.

"""
struct LogClosureAmplitudes{K} <: VLBIDataProducts{K}
    keywords::K
end

"""
    Visibilities(;kwargs...)

Type to specify to extract the log closure amplitudes table in the [`extract_table`](@ref) function.
Optional keywords are passed through `extract_table` to specify additional option.

# Special keywords for eht-imaging with Pyehtim.jl
For a list of potential keyword arguments see `eht-imaging` and `add_amp` command for obsdata.
"""
struct VisibilityAmplitudes{K} <: VLBIDataProducts{K}
    keywords::K
end

"""
    Visibilities(;kwargs...)

Type to specify to extract the complex visibilities table in the [`extract_table`](@ref) function.
Optional keywords are passed through `extract_table` to specify additional option.

# Special keywords for eht-imaging with Pyehtim.jl
Any keyword arguments are ignored for now. Use eht-imaging directly to modify the data.
"""
struct Visibilities{K} <: VLBIDataProducts{K}
    keywords::K
end

"""
    Coherencies(;kwargs...)

Type to specify to extract the coherency matrices table in the [`extract_table`](@ref) function.
Optional keywords are passed through `extract_table` to specify additional option.

# Special keywords for eht-imaging with Pyehtim.jl
Any keyword arguments are ignored for now. Use eht-imaging directly to modify the data.
"""
struct Coherencies{K} <: VLBIDataProducts{K}
    keywords::K
end

for c in [:ClosurePhases, :LogClosureAmplitudes, :VisibilityAmplitudes, :Visibilities, :Coherencies]
    @eval begin
        $(c)(; kwargs...) = $(c)(kwargs)
    end
end

"""
    extract_table(obs, dataproducts::VLBIDataProducts)

Extract an [`Comrade.EHTObservationTable`](@ref) table of data products `dataproducts`.
To pass additional keyword for the data products you can pass them as keyword arguments
to the data product type. For a list of potential data products see `subtypes(Comrade.VLBIDataProducts)`.

# Example
```julia-repl
julia> dlcamp, dcphase = extract_table(obs, LogClosureAmplitudes(;snrcut=3.0), ClosurePhases(;snrcut=3.0, cut_trivial=true))
julia> dcoh = extract_table(obs, Coherencies())
julia> dvis = extract_table(obs, VisibilityAmplitudes())
```
"""
function extract_table(obs, dataproducts::VLBIDataProducts...)
    @assert length(dataproducts) >= 1 "No dataproducts passed to `extract_table`"
    return map(x -> extract_table(obs, x), dataproducts)
end

function extract_table(obs, dataproduct::ClosurePhases)
    return extract_cphase(obs; keywords(dataproduct)...)
end

function extract_table(obs, dataproduct::LogClosureAmplitudes)
    return extract_lcamp(obs; keywords(dataproduct)...)
end

function extract_table(obs, dataproduct::Visibilities)
    return extract_vis(obs; keywords(dataproduct)...)
end

function extract_table(obs, dataproduct::VisibilityAmplitudes)
    return extract_amp(obs; keywords(dataproduct)...)
end

function extract_table(obs, dataproduct::Coherencies)
    return extract_coherency(obs; keywords(dataproduct)...)
end

# internal methods to extract information from `obs`
"""
    extract_cphase(obs; kwargs...)

Extracts the closure phases from an `obs`.
This is an internal method for dispatch. Only use this if
interfacing Comrade with a new data type.
"""
function extract_cphase    end
"""
    extract_lcamp(obs; kwargs...)

Extracts the log-closure amplitudes from an `obs`.
This is an internal method for dispatch. Only use this if
interfacing Comrade with a new data type.
"""
function extract_lcamp     end
"""
    extract_amp(obs; kwargs...)

Extracts the visibility amplitudes from an `obs`.
This is an internal method for dispatch. Only use this if
interfacing Comrade with a new data type.
"""
function extract_amp       end
"""
    extract_vis(obs; kwargs...)

Extracts the stokes I complex visibilities from an obs.
This is an internal method for dispatch. Only use this if
interfacing Comrade with a new data type.
"""
function extract_vis       end
"""
    extract_coherency(obs; kwargs...)

Extracts the full coherency matrix from an observation.
This is an internal method for dispatch. Only use this if
interfacing Comrade with a new data type.
"""
function extract_coherency end
