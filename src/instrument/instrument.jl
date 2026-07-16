using LinearAlgebra
using SparseArrays
import Distributions
using Statistics
using PrettyTables

"""
    Segment(lo, hi)

A half-open axis segment `[lo, hi)` (time in UTC hours, or frequency in Hz): the
construction-time currency of the segmentation layer ([`timestamps`](@ref) and
`freqchannels` return grids of these). Only the explicit edges are stored — the center
is derived — so datum→segment membership never depends on reconstructing edges from
center/width arithmetic. Segments are transient: the parameter store keeps their centers
and edges as plain axis vectors (see `SiteLookup`).
"""
struct Segment{T <: Real}
    lo::T
    hi::T
end

_center(s::Segment) = (s.lo + s.hi) / 2
_width(s::Segment) = s.hi - s.lo
Base.in(x::Number, s::Segment) = s.lo ≤ x < s.hi
Base.isless(a::Segment, b::Segment) = a.lo < b.lo
Base.isless(x::Number, s::Segment) = x < s.lo
# a segment sorts before `x` iff it lies entirely below it (half-open, so `hi ≤ x`);
# `searchsortedfirst(segments, x)` then lands on the first segment that can contain `x`
Base.isless(s::Segment, x::Number) = s.hi ≤ x
Base.Broadcast.broadcastable(s::Segment) = Ref(s)


include("site_array.jl")
include("feedrotations.jl")
include("jonesmatrices.jl")
include("priors/priors.jl")
include("instrument_transforms.jl")
include("model.jl")
include("macro.jl")
include("utility.jl")
include("caltable.jl")
