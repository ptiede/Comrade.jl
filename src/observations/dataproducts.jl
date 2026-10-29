export ClosurePhases, LogClosureAmplitudes


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

"""
    ClosurePhases()

The kind of a [`ClosureTable`](@ref) of closure phases; see [`closure_phases`](@ref).
"""
struct ClosurePhases end

"""
    LogClosureAmplitudes()

The kind of a [`ClosureTable`](@ref) of log closure amplitudes; see
[`log_closure_amplitudes`](@ref).
"""
struct LogClosureAmplitudes end

"""
    load_array_txt(path)

Read an ehtim-style antenna text file and return a dictionary keyed by site
symbol. Each value is a NamedTuple with the columns present in the file:
`(SEFD1, SEFD2, fr_parallactic, fr_elevation, fr_offset)` plus the d-term
fields `(DR, DL)` if present. `fr_offset` is converted from degrees (the
ehtim file convention) to radians.

The file must begin with a commented header naming the columns, e.g.

    # site  X  Y  Z  SEFDR  SEFDL  DR_RE  DR_IM  DL_RE  DL_IM  FR_PAR  FR_ELEV  FR_OFF

Columns are matched by name (case-insensitive), so the order of fields and the
relative position of the D-term and FR blocks are flexible. `X, Y, Z` are
ignored. The result can be passed to [`reset_mounts!`](@ref).
"""
function load_array_txt(path::AbstractString)
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
