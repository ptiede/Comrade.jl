export extract_table,
    ClosurePhases, LogClosureAmplitudes,
    VisibilityAmplitudes, Visibilities, Coherencies


# Traits that decide the domain of the data. We include this to prevent additional computation if
# e.g., on visibility data is considered.
abstract type ComradeDataType end
struct NoData <: ComradeDataType end
struct DualData <: ComradeDataType end
struct VisData <: ComradeDataType end
struct ImgData <: ComradeDataType end

datatype(::Type, ::Type{<:Nothing}) = VisData()
datatype(::Type{<:Nothing}, ::Type) = ImgData()
datatype(::Type{<:Nothing}, ::Type{<:Nothing}) = throw(ArgumentError("No data in uv plane or image plane is provided"))
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

"""
    load_array_txt(path) -> Dict{Symbol, NamedTuple}

Read an ehtim-style antenna text file and return a dictionary keyed by site
symbol. Each value is a NamedTuple with the columns present in the file:
`(SEFD1, SEFD2, fr_parallactic, fr_elevation, fr_offset)` plus the d-term
fields `(DR, DL)` if present. `fr_offset` is converted from degrees (the
ehtim file convention) to radians.

The file must begin with a commented header naming the columns, e.g.

    # site  X  Y  Z  SEFDR  SEFDL  DR_RE  DR_IM  DL_RE  DL_IM  FR_PAR  FR_ELEV  FR_OFF

Columns are matched by name (case-insensitive), so the order of fields and the
relative position of the D-term and FR blocks are flexible. `X, Y, Z` are
ignored here — antenna positions come from the uvfits AIPS-AN table. The
result can be passed to any `extract_*` as
`array_overrides=load_array_txt("array.txt")`.
"""
function Comrade.load_array_txt(path::AbstractString)
    # Aliases for column names we care about. Matching is case-insensitive.
    aliases = Dict(
        :site => ("site", "name", "station"),
        :sefdr => ("sefdr", "sefd1", "sefd_r"),
        :sefdl => ("sefdl", "sefd2", "sefd_l"),
        :dr_re => ("dr_re", "dre", "drre", "d_r_re"),
        :dr_im => ("dr_im", "dim", "drim", "d_r_im"),
        :dl_re => ("dl_re", "dlre", "d_l_re"),
        :dl_im => ("dl_im", "dlim", "d_l_im"),
        :fr_par => ("fr_par", "frpar", "fr_parallactic", "f_par"),
        :fr_el => ("fr_elev", "frelev", "fr_elevation", "f_el", "fr_el"),
        :fr_off => ("fr_off", "fr_offset", "froff", "f_off"),
    )

    header = nothing
    rows = String[]
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "#")
            # Use the *last* commented line before data as the header.
            header = strip(lstrip(line, '#'))
        else
            push!(rows, line)
        end
    end
    header === nothing && error(
        "load_array_txt: $path has no '#'-commented header line — column layout cannot be inferred."
    )

    header_toks = split(header)
    # idx[:key] = column index in the data rows
    idx = Dict{Symbol, Int}()
    for (key, names) in aliases
        for (i, tok) in enumerate(header_toks)
            if lowercase(tok) in names
                idx[key] = i
                break
            end
        end
    end

    required = (:site, :sefdr, :sefdl, :fr_par, :fr_el, :fr_off)
    missing_cols = filter(k -> !haskey(idx, k), collect(required))
    isempty(missing_cols) || error(
        "load_array_txt: header is missing required columns $missing_cols. Header was: $header"
    )

    has_dterms = all(haskey(idx, k) for k in (:dr_re, :dr_im, :dl_re, :dl_im))
    ncols_needed = maximum(values(idx))

    overrides = Dict{Symbol, NamedTuple}()
    for line in rows
        toks = split(line)
        length(toks) >= ncols_needed || continue
        site = Symbol(toks[idx[:site]])
        sefd_r = parse(Float64, toks[idx[:sefdr]])
        sefd_l = parse(Float64, toks[idx[:sefdl]])
        fr_par = parse(Float64, toks[idx[:fr_par]])
        fr_el = parse(Float64, toks[idx[:fr_el]])
        fr_off = deg2rad(parse(Float64, toks[idx[:fr_off]]))
        nt = (
            SEFD1 = sefd_r, SEFD2 = sefd_l,
            fr_parallactic = fr_par, fr_elevation = fr_el, fr_offset = fr_off,
        )
        if has_dterms
            dr = complex(
                parse(Float64, toks[idx[:dr_re]]),
                parse(Float64, toks[idx[:dr_im]]),
            )
            dl = complex(
                parse(Float64, toks[idx[:dl_re]]),
                parse(Float64, toks[idx[:dl_im]]),
            )
            nt = merge(nt, (DR = dr, DL = dl))
        end
        overrides[site] = nt
    end
    return overrides
end


