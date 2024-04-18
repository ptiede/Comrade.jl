export ScanSeg, IntegSeg, TrackSeg

"""
    $(TYPEDEF)

The data segmentation scheme to use. This is used specify how often we want various
instrument hyperparameters to change. A user subtyping this expression must implement
the following functions:

 - [`timestamps`](@ref): Specifies the time region for each segmentation scheme given an array
"""
abstract type ObsSegmentation end

# Track is for quantities that remain static across an entire observation
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over a `track`, i.e., the observation "night".
"""
struct TrackSeg <: ObsSegmentation end

# Scan is for quantities that are constant across a scan
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over a `scan`.
"""
struct ScanSeg <: ObsSegmentation end


# Integration is for quantities that change every integration time
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over the time stamps in the data.
If the data is scan-averaged before then `IntegSeg` will be identical to `ScanSeg`.
"""
struct IntegSeg <: ObsSegmentation end

"""
    timestamps(seg::Segmentation, array::AbstractArrayConfiguration)

Return the time stamps or really a vector of integration time regions for a
given segmentation scheme `seg` and array configuration `array`.
"""
function timestamps end

function timestamps(::ScanSeg, array)
    st     = array.scans
    scanid = 1:length(st)
    mjd    = array.mjd

    # Shift the central time to the middle of the scan
    dt = (st.stop .- st.start)
    t0 = st.start .+ dt./2
    return IntegrationTime.(mjd, scanid, t0, dt)
end

getscan(scans, t) = findfirst(i->scans.start[i]≤t<scans.stop[i], 1:length(scans))


function timestamps(::IntegSeg, array)
    ts = unique(array[:T])
    st     = array.scans
    scanids = getscan.(Ref(st), ts)
    mjd    = array.mjd

    # TODO build in the dt into the data format
    dt = minimum(diff(ts))
    return IntegrationTime.(mjd, scanids, ts, dt)
end

function timestamps(::TrackSeg, array)
    st     = array.scans
    mjd    = array.mjd

    tstart, tend = extrema(array[:T])
    dt = tend - tstart

    # TODO build in the dt into the data format
    return (IntegrationTime(mjd, 1:length(st), (tend-tstart)/2 + tstart, dt),)
end
