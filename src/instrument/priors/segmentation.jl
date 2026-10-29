export ScanSeg, IntegSeg, TrackSeg, EpochSeg

"""
    $(TYPEDEF)

The data segmentation scheme to use. This is used specify how often we want various
instrument hyperparameters to change. A user subtyping this expression must implement
the following functions:

 - [`timestamps`](@ref): Specifies the time region for each segmentation scheme given a visibility table
"""
abstract type Segmentation end

abstract type TimeSegmentation <: Segmentation end

# Track is for quantities that remain static across an entire observation
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over the whole table: every track of a
table joined from several nights shares one value. See [`EpochSeg`](@ref) for one value per
track.
"""
struct TrackSeg <: TimeSegmentation end

"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over each track of `info(t).tracks`, an
observing session such as one night, and changes from track to track.
"""
struct EpochSeg <: TimeSegmentation end

# Scan is for quantities that are constant across a scan
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over a `scan`.
"""
struct ScanSeg <: TimeSegmentation end


# Integration is for quantities that change every integration time
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over the time stamps in the data.
If the data is scan-averaged before then `IntegSeg` will be identical to `ScanSeg`.
"""
struct IntegSeg <: TimeSegmentation end

"""
    timestamps(seg::TimeSegmentation, t::VisTable)

The time intervals, as `IntegrationTime`s, into which the segmentation `seg` divides the
measured values of the visibility table `t`.
"""
function timestamps end

function timestamps(::ScanSeg, t::VisTable)
    st = info(t).scans
    dt = (st.stop .- st.start)
    t0 = st.start .+ dt ./ 2
    return IntegrationTime.(_mjd(t), t0, dt)
end

function timestamps(::IntegSeg, t::VisTable)
    ts = sort(unique(_measuredpoints(t).Ti))
    # TODO build in the dt into the data format
    dt = length(ts) <= 1 ? 0.1 / 3600 : minimum(diff(ts))
    return IntegrationTime.(_mjd(t), ts, dt)
end

function timestamps(::TrackSeg, t::VisTable)
    tstart, tend = extrema(_measuredpoints(t).Ti)
    tend = tend + 1
    dt = 2 * (tend - tstart)
    if iszero(dt)
        dt = 2 * 1 / 3600
    end
    return (IntegrationTime(_mjd(t), (tend - tstart) / 2 + tstart, dt),)
end

# One segment per track, centered on it and widened by an hour, or by half the gap to the
# nearest track if that is shorter, so it holds the time stamps at the track's ends.
function timestamps(::EpochSeg, t::VisTable)
    tr = info(t).tracks
    gap = minimum(tr.start[2:end] .- tr.stop[1:(end - 1)]; init = Inf)
    pad = min(1.0, gap / 2)
    return map((a, b) -> IntegrationTime(_mjd(t), (a + b) / 2, b - a + 2 * pad), tr.start, tr.stop)
end

# The frequency segment of every instrument parameter: one spanning all measured channels of
# `t`, so a station's parameter is shared by its channels. Frequency dependence belongs in the
# Jones function, which reads each cell's frequency `c.Fr`.
function bandchannel(t::VisTable)
    ch = info(t).channels
    lo, hi = extrema(_measuredpoints(t).Fr)
    width(f) = ch.width[findfirst(==(f), ch.frequency)]
    a, b = lo - width(lo) / 2, hi + width(hi) / 2
    return FrequencyChannel((a + b) / 2, b - a, 1)
end
