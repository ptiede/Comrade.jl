"""
    $(TYPEDEF)
Abstract trait for Fourier transform trait. This will decide whether
a given model should use the analytic (if it exists) or numerical
FT.
"""
abstract type VisStyle end

struct VisAnalytic <: VisStyle end
struct VisNumeric <: VisStyle end

"""
    VisStyle(::AbstractModel)
Sets the trait for the Fourier transform to be used for a given model.
The default value is IsNumeric() which says to use the numerical Fourier
transform.
"""
VisStyle(::Type{<:AbstractModel}) = VisNumeric()






abstract type ImStyle end

struct ImAnalytic <: ImStyle end
struct ImNumeric <: ImStyle end

"""
    ImStyle(::Type{AbstractModel})
Determines whether the image is pointwise analytic, i.e. we can evaluate
its intensity as a pixel location.

If the model has the IsAnalytic() trait then the user must implement a function
intensity, that computes the intensity at a location.

If the model has the IsNumeric() trait then the user must specify intensitymap! that
return the entire image on a grid. The intensity function then works by memoizing this function
and interpolating on it.
"""
ImStyle(::Type{<:AbstractModel}) = ImAnalytic()
