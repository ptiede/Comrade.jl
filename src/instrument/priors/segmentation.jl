export ScanSeg, IntegSeg, TrackSeg, SpectralWindow, FullBand

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

Return the grid of half-open time [`Segment`](@ref Comrade.Segment)s implied by the
segmentation scheme `seg` and array configuration `array`.
"""
function timestamps end

timestamps(::ScanSeg, array) = Segment.(array.scans.start, array.scans.stop)

function timestamps(::IntegSeg, array)
    ts = sort(unique(array[:Ti]))

    # TODO build in the dt into the data format
    if length(ts) <= 1
        # arbritrarily set the dt to 1
        dt = 0.1 / 3600
    else
        dt = minimum(diff(ts))
    end
    return Segment.(ts .- dt / 2, ts .+ dt / 2)
end

function timestamps(::TrackSeg, array)
    # A track must contain every finer segmentation, and a scan brackets its (averaged)
    # datum times — so span the scan table as well as the data, otherwise a ScanSeg
    # segment (and its center) can start before the track does and segment-membership
    # matching across the two segmentations fails.
    tstart, tend = float.(extrema(array[:Ti]))
    scans = array.scans
    if !isempty(scans)
        tstart = min(tstart, float(minimum(scans.start)))
        tend = max(tend, float(maximum(scans.stop)))
    end
    return (Segment(tstart, nextfloat(tend)),)
end


"""
    $(TYPEDEF)

Frequency segmentation with one independent parameter per observed frequency channel.
This is the default `freqseg` of the site priors.
"""
struct SpectralWindow <: FrequencySegmentation end


function freqchannels(::SpectralWindow, array)
    Fr = float.(sort(unique(array[:Fr])))
    bw = float(array.bandwidth)
    # a channel is centered on its observed frequency; guard the half-open upper edge for
    # degenerate (zero) bandwidths
    return Segment.(Fr .- bw / 2, max.(Fr .+ bw / 2, nextfloat.(Fr)))
end

"""
    $(TYPEDEF)

Frequency segmentation with a single parameter shared across the whole observed band:
one channel covering every observed frequency. Use as the `freqseg` of a site prior for
quantities that do not vary across IFs, e.g.
`IIDSitePrior(ScanSeg(), dist; freqseg = FullBand())` gives one gain per scan shared by
all frequency channels. A frequency-dependent `param_map` still sees each *datum's*
frequency, so spectral gain laws compose with a `FullBand` gain. Reference-antenna
selection treats a band-spanning parameter as a candidate in every per-channel reference
group it overlaps.
"""
struct FullBand <: FrequencySegmentation end

function freqchannels(::FullBand, array)
    lo, hi = float.(extrema(array[:Fr]))
    return (Segment(lo, nextfloat(hi)),)
end
