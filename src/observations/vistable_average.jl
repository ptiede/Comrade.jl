export average, ByScan, ChannelBins, ByChannel, averaging_axis, times

"""
    times(t::VisTable)

The distinct time stamps of `t`, in hours since `info(t).epoch`, sorted: the `Ts` lookup of
its station grid.
"""
times(t::VisTable) = collect(DD.lookup(stations(t), Ts))

"""
    average(t::VisTable, by...)

`t` with the values sharing a label collapsed into one, along the axis each criterion `by`
labels ([`averaging_axis`](@ref)): time stamps for a time criterion, channels for a channel
criterion ([`ChannelBins`](@ref), [`ByChannel`](@ref)). `by(t)` gives one label per time stamp
([`times`](@ref)) or per channel, or a single label for the whole axis. Several criteria apply
in turn, the time criteria first; the order does not change the result.

```julia
average(t, ByScan())                                  # one record per baseline and scan
average(t, t -> floor.(Int, times(t) ./ (10 / 60)))   # one per baseline and ten minutes
average(t, ByScan(), ChannelBins(4))                  # per scan, four channels to one
```

Each averaged value is the inverse-variance weighted mean `Σ wV / Σ w` of the measured values
under its label and its weight their sum `Σ w`; where none is measured it is unmeasured.

Over time, the records of each baseline under a label become one record. Its time is the
mid-point of the label's first and last time stamps and its `u`, `v` the means of its
records'; neither depends on the weights. The station grid is recomputed from the station
positions at the new time stamps.

Over channels, an averaged channel's frequency is the mean of its channels' frequencies and
its width their summed widths. The channels under a label need not be adjacent; the averaged
channels are sorted by frequency.
"""
function average(t::VisTable, by, more...)
    eltype(measurement(t)) <: Real &&
        throw(ArgumentError("average the visibilities, then take their amplitudes"))
    criteria = (by, more...)
    order = sortperm(collect(criteria); by = c -> averaging_axis(c) isa Ti ? 0 : 1)
    return foldl((x, i) -> _average(x, criteria[i], averaging_axis(criteria[i])), order; init = t)
end

"""
    averaging_axis(by)

The axis the criterion `by` labels for [`average`](@ref): `Ti()`, the default, so that any
function of a table is a time criterion, or `Fr()`. A criterion type of your own that labels
channels extends it.
"""
averaging_axis(by) = Ti()

"""
    ByScan()

The time criterion for [`average`](@ref) that labels each time stamp with the scan of
`info(t).scans` holding it. A time stamp outside every scan is an error.
"""
struct ByScan end

function (::ByScan)(t::VisTable)
    sc = info(t).scans
    return map(times(t)) do τ
        s = findfirst(i -> sc.start[i] <= τ <= sc.stop[i], eachindex(sc.start))
        isnothing(s) && throw(ArgumentError("the time stamp $τ h is outside every scan"))
        s
    end
end

"""
    ChannelBins(n)

The channel criterion for [`average`](@ref) that gathers each run of `n` channels, in order of
frequency, into one. The number of channels must be a multiple of `n`.
"""
struct ChannelBins
    n::Int
    function ChannelBins(n::Integer)
        n > 0 || throw(ArgumentError("`ChannelBins` gathers a positive number of channels, got $n"))
        return new(n)
    end
end

averaging_axis(::ChannelBins) = Fr()

function (b::ChannelBins)(t::VisTable)
    nch = length(info(t).channels)
    nch % b.n == 0 ||
        throw(ArgumentError("cannot gather $nch channels in runs of $(b.n); the run length must divide the number of channels"))
    return (0:(nch - 1)) .÷ b.n
end

"""
    ByChannel(f)

A channel criterion for [`average`](@ref) from the function `f`: `f(t)` gives one label per
channel of `t`, or a single label for all of them.

```julia
average(t, ByChannel(t -> round.(Int, info(t).channels.frequency ./ 8.0e6)))   # 8 MHz bins
```
"""
struct ByChannel{F}
    f::F
end

averaging_axis(::ByChannel) = Fr()
(b::ByChannel)(t::VisTable) = b.f(t)

# One label per element of an axis of length `n`, from a criterion's labels.
function _labels(given, n, what)
    given isa AbstractArray || return fill(given, n)
    length(given) == n ||
        throw(DimensionMismatch("a criterion gave $(length(given)) labels for $n $what"))
    return collect(given)
end

# The members of each label, in the order the labels first occur.
function _members(labels)
    order = unique(labels)
    index = Dict(l => i for (i, l) in pairs(order))
    members = [Int[] for _ in order]
    for (i, l) in pairs(labels)
        push!(members[index[l]], i)
    end
    return members, [index[l] for l in labels]
end

_safediv(x, w) = x ./ ifelse.(w .> 0, w, one(eltype(w)))

function _average(t::VisTable, by, ::Ti)
    d = domain(t)
    τ = times(t)
    tlabel = _labels(by(t), length(τ), "time stamps")
    tmembers, tgroup = _members(tlabel)
    keys = [(tgroup[d.ts[n]], d.ant1[n], d.ant2[n]) for n in eachindex(d.ts)]
    members, group = _members(keys)
    n, ng = length(keys), length(members)
    A = sparse(group, 1:n, 1.0, ng, n)

    V, w = _dense(measurement(t)), _dense(weight(t))
    rest = size(V)[2:end]
    W = reshape(w, n, :)
    Ws = A * W
    Vs = _safediv(A * (W .* reshape(V, n, :)), Ws)

    c = ComradeBase.coords(d)
    nrecords = A * ones(n)
    mid = [(τ[first(m)] + τ[last(m)]) / 2 for m in tmembers]
    rep = first.(members)
    return VisTable(
        t.pol, reshape(Vs, ng, rest...); weight = reshape(Ws, ng, rest...),
        u = (A * c.u) ./ nrecords, v = (A * c.v) ./ nrecords, Ti = mid[first.(keys[rep])],
        ant1 = d.ant1[rep], ant2 = d.ant2[rep], info = info(t),
    )
end

function _average(t::VisTable, by, ::Fr)
    ch = info(t).channels
    members, _ = _members(_labels(by(t), length(ch), "channels"))
    frequency = [sum(ch.frequency[m]) / length(m) for m in members]
    order = sortperm(frequency)
    members = members[order]

    V, w = _dense(measurement(t)), _dense(weight(t))
    gather(x, m) = sum(selectdim(x, 2, m); dims = 2)
    Ws = cat((gather(w, m) for m in members)...; dims = 2)
    Vs = _safediv(cat((gather(w .* V, m) for m in members)...; dims = 2), Ws)

    a = info(t)
    channels = (; frequency = frequency[order], width = [sum(ch.width[m]) for m in members])
    newinfo = ArrayInfo(; stations = a.stations, channels, scans = a.scans, tracks = a.tracks, source = a.source, epoch = a.epoch)
    d = domain(t)
    c = ComradeBase.coords(d)
    g = stations(t)
    return VisTable(
        t.pol, Vs; weight = Ws, u = c.u, v = c.v, Ti = d.Ti, ant1 = d.ant1, ant2 = d.ant2,
        info = newinfo, elevation = parent(g.elevation), parallactic = parent(g.parallactic),
    )
end
