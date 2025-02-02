export flag, select_baseline, add_fractional_noise!, add_fractional_noise

function add_fractional_noise!(dvis, ferr)
    map!(dvis[:noise], dvis[:noise], dvis[:measurement]) do e, m
        fe =  sqrt.(e.^2 .+ ferr.^2*abs2(tr(m))/2)
        return fe
    end
    return dvis
end

function add_fractional_noise(dvis, ferr)
    dvis = deepcopy(dvis)
    return add_fractional_noise!(dvis, ferr)
end

"""
    flag(condition, data)

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
function flag(condition, dvis; negate=false)
    inds = findall(x->condition(x), datatable(dvis))
    return dvis[inds]
end

function Base.filter(f, dvis::EHTObservationTable)
    return flag(x->!f(x), dvis)
end



function select_baseline(dvis, bl::Tuple{Symbol, Symbol})
    cond(x) = !((x.baseline.sites == bl) || (x.baseline.sites == reverse(bl)))
    return filter(cond, dvis) 
end


"""
    residual_data(vis::AbstractArray, data::EHTObservationTable)

Compute the residuals for the model visibilities `vis` and the data `data`.
The residuals are not normalized and the returned object is an `EHTObservationTable`.
"""
function residual_data(vis, data::EHTObservationTable{A}) where {A<:EHTClosurePhaseDatum}
    phase = measurement(data)
    err = noise(data)

    mphase = closure_phases(vis, designmat(arrayconfig(data)))
    res = @. atan(sin(phase - mphase), cos(phase - mphase))
    return EHTObservationTable{A}(res, err, arrayconfig(data))
end


function residual_data(vis, dlca::EHTObservationTable{A}) where {A<:EHTLogClosureAmplitudeDatum}
    phase = measurement(dlca)
    err = noise(dlca)
    mphase = logclosure_amplitudes(vis, designmat(arrayconfig(dlca)))
    res = (phase .- mphase)
    return EHTObservationTable{A}(res, err, arrayconfig(dlca))
end

function residual_data(vis, damp::EHTObservationTable{A}) where {A<:EHTVisibilityAmplitudeDatum}
    mamp = abs.(vis)
    amp = measurement(damp)
    res = (amp - mamp)
    return EHTObservationTable{A}(res, noise(damp), arrayconfig(damp))
end

function residual_data(vis, dvis::EHTObservationTable{A}) where {A<:EHTCoherencyDatum}
    coh = measurement(dvis)
    res = coh .- vis
    return EHTObservationTable{A}(res, noise(dvis), arrayconfig(dvis))
end

function residual_data(mvis, dvis::EHTObservationTable{A}) where {A<:EHTVisibilityDatum}
    vis = measurement(dvis)
    res = (vis - mvis)
    return EHTObservationTable{A}(res, noise(dvis), arrayconfig(dvis))
end
