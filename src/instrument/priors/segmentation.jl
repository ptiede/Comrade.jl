export ScanSeg, IntegSeg, TrackSeg

"""
    $(TYPEDEF)

The data segmentation scheme to use. This is used specify how often we want various
instrument hyperparameters to change. A user subtyping this expression must implement
the following functions:

 - [`timestamps`](@ref): Specifies the time region for each segmentation scheme given an array
"""
abstract type Segmentation end

abstract type TimeSegmentation <: Segmentation end

abstract type FrequencySegmentation <: Segmentation end

# Track is for quantities that remain static across an entire observation
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over a `track`, i.e., the observation "night".
"""
struct TrackSeg <: TimeSegmentation end

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
    timestamps(seg::TimeSegmentation, array::AbstractArrayConfiguration)

Return the time stamps or really a vector of integration time regions for a
given segmentation scheme `seg` and array configuration `array`.
"""
function timestamps end

function timestamps(::ScanSeg, array)
    st = array.scans
    mjd = array.mjd
    # Shift the central time to the middle of the scan
    dt = (st.stop .- st.start)
    t0 = st.start .+ dt ./ 2

    return IntegrationTime.(mjd, t0, dt)
end

# getscan(scans, t) = findfirst(i->scans.start[i]≤t<scans.stop[i], 1:length(scans))


function timestamps(::IntegSeg, array)
    ts = unique(array[:Ti])
    mjd = array.mjd

    # TODO build in the dt into the data format
    if length(ts) <= 1
        # arbritrarily set the dt to 1
        dt = 1 / 3600
    else
        dt = minimum(diff(ts))
    end
    return IntegrationTime.(mjd, ts, dt)
end

function timestamps(::TrackSeg, array)
    mjd = array.mjd

    tstart, tend = extrema(array[:Ti])
    tend = tend + 1
    dt = tend - tstart
    if iszero(dt)
        dt = 1 / 3600
    end
    # TODO build in the dt into the data format
    return (IntegrationTime(mjd, (tend - tstart) / 2 + tstart, dt),)
end


struct SpectralWindow <: FrequencySegmentation end


function freqchannels(::SpectralWindow, array)
    Fr = unique(array[:Fr])
    if length(Fr) > 1
        Fr[end] = nextfloat(Fr[end])
    end
    return FrequencyChannel.(Fr, array.bandwidth, 1:length(Fr))
end
