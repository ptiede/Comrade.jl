export ScanSeg, IntegSeg, TrackSeg

"""
    $(TYPEDEF)

The data segmentation scheme to use. This is important for constructing a [`JonesCache`](@ref)
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

## Warning
Currently we do not explicity track the telescope scans. This will be fixed in a future version.
Right now `ScanSeg` and `TrackSeg` are the same
"""
struct ScanSeg <: ObsSegmentation end


# Integration is for quantities that change every integration time
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over a correlation integration.
"""
struct IntegSeg <: ObsSegmentation end
