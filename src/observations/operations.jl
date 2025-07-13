export flag, select_baseline, add_fractional_noise!, add_fractional_noise

function add_fractional_noise!(dvis::EHTObservationTable{<:EHTCoherencyDatum}, ferr)
    map!(dvis[:noise], dvis[:noise], dvis[:measurement]) do e, m
        # we do it like this so I don't NaN everything
        # if a feed is missing
        if isnan(m[1, 1])
            err = m[2, 2]
        elseif isnan(m[2, 2])
            err = m[1, 1]
        else
            err = tr(m) / 2
        end
        fe = sqrt.(e .^ 2 .+ ferr .^ 2 * abs2(err))
        return fe
    end
    return dvis
end

function add_fractional_noise!(
        dvis::EHTObservationTable{<:Union{EHTVisibilityDatum, EHTVisibilityAmplitudeDatum}},
        ferr
    )
    map!(dvis[:noise], dvis[:noise], dvis[:measurement]) do e, m
        fe = sqrt.(e .^ 2 .+ ferr .^ 2 * abs2(m))
        return fe
    end
    return dvis
end

function add_fractional_noise!(dcl::EHTObservationTable{<:ClosureProducts}, ferr)
    conf = arrayconfig(dcl)
    nois = noise(dcl)
    vis = getfield(conf, :vis)
    nvi = getfield(conf, :noise)
    dmat = designmat(conf)
    map!(nvi, nvi, vis) do e, m
        fe = sqrt.(e .^ 2 .+ ferr .^ 2 * abs2(m))
        return fe
    end
    # update the noise covariance matrix
    nois .= dmat * Diagonal(abs2.(nvi ./ vis)) * transpose(dmat)
    return dcl
end


"""
    add_fractional_noise(dvis::AbstractObservationTable, ferr)

Adds fractional noise to the data in `dvis` with the fractional error `ferr`, i.e. the returned
data is identical except the noise is replaced with

   √(σ² + ferr² * |V|²)

where `σ` is dvis[:noise] and `V` is dvis[:measurement].

!!! warning
    When the datum type <: ClosureProducts, we do not update the noise covariance matrix.
    This means that cross-correlation may be slightly incorrect.
"""
function add_fractional_noise(dvis, ferr)
    dvis = deepcopy(dvis)
    return add_fractional_noise!(dvis, ferr)
end

"""
    flag(condition, data::AbstractObservationTable)

Flags the data in `vis` that satisfy the condition `condition`.
For the opposite effect call `filter(condition, data)`

`condition` is a function that accepts a tuple (baseline, measurement, noise) and returns a boolean.
`baseline` are the fields that define the measurement like in `datatable(vis).baseline`.

# Example
To flag all data with a uvdistance less than 0.1e9 use:
```julia
cond(x) = uvdist(x) < 0.1e9
flag_data(vis, cond)
```
"""
function flag(condition, dvis::AbstractObservationTable)
    inds = findall(x -> !condition(x), datatable(dvis))
    return dvis[inds]
end

"""
    filter(condition, data::AbstractObservationTable)

Filters all the data in `data` that satisfy the condition `condition`, i.e.
if condition returns true for a given datum, it is kept in the returned object.

This is the converse of `flag`.
"""
function Base.filter(f, dvis::AbstractObservationTable)
    return flag(x -> !f(x), dvis)
end


"""
    select_baseline(data::AbstractObservationTable, bl::Tuple{Symbol, Symbol})

Selects the baseline `bl` from the `data` object. The baseline is a tuple of symbols,
but the order of the symbols does not matter.
"""
function select_baseline(dvis, bl::Tuple{Symbol, Symbol})
    cond(x) = ((x.baseline.sites == bl) || (x.baseline.sites == reverse(bl)))
    fl = filter(cond, dvis)
    dt = map(datatable(fl)) do d
        si = d.baseline.sites
        si != bl && return flipbaseline(d)
        return d
    end
    return rebuild(fl, dt)
end

function flipbaseline(datum::EHTCoherencyDatum)
    return EHTCoherencyDatum(
        datum.measurement',
        datum.noise',
        flipbaseline(datum.baseline)
    )
end

function flipbaseline(datum::EHTVisibilityDatum)
    return EHTCoherencyDatum(
        conj(datum.measurement),
        conj(datum.noise),
        flipbaseline(datum.baseline)
    )
end

function flipbaseline(datum::EHTVisibilityAmplitudeDatum)
    return datum
end


"""
    residual_data(vis::AbstractArray, data::EHTObservationTable)

Compute the residuals for the model visibilities `vis` and the data `data`.
The residuals are not normalized and the returned object is an `EHTObservationTable`.
"""
function residual_data(vis, data::EHTObservationTable{A}) where {A <: EHTClosurePhaseDatum}
    phase = measurement(data)
    err = noise(data)

    mphase = closure_phases(vis, designmat(arrayconfig(data)))
    res = @. atan(sin(phase - mphase), cos(phase - mphase))
    return EHTObservationTable{A}(res, err, arrayconfig(data))
end


function residual_data(vis, dlca::EHTObservationTable{A}) where {A <: EHTLogClosureAmplitudeDatum}
    phase = measurement(dlca)
    err = noise(dlca)
    mphase = logclosure_amplitudes(vis, designmat(arrayconfig(dlca)))
    res = (phase .- mphase)
    return EHTObservationTable{A}(res, err, arrayconfig(dlca))
end

function residual_data(vis, damp::EHTObservationTable{A}) where {A <: EHTVisibilityAmplitudeDatum}
    mamp = abs.(vis)
    amp = measurement(damp)
    res = (amp - mamp)
    return EHTObservationTable{A}(res, noise(damp), arrayconfig(damp))
end

function residual_data(vis, dvis::EHTObservationTable{A}) where {A <: EHTCoherencyDatum}
    coh = measurement(dvis)
    res = coh .- vis
    return EHTObservationTable{A}(res, noise(dvis), arrayconfig(dvis))
end

function residual_data(mvis, dvis::EHTObservationTable{A}) where {A <: EHTVisibilityDatum}
    vis = measurement(dvis)
    res = (vis - mvis)
    return EHTObservationTable{A}(res, noise(dvis), arrayconfig(dvis))
end
