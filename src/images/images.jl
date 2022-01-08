abstract type AbstractIntensityMap{T,S} <: AbstractMatrix{T} end

abstract type AbstractPolarizedMap{I,Q,U,V} end


"""
Pulse
Pixel response function for a radio image model. This makes
a discrete sampling continuous by picking a certain *smoothing*
kernel for the image.

# Notes
To see the implemented Pulses please use the subtypes function i.e.
`subtypes(Pulse)`
"""
abstract type Pulse <: AbstractModel end

visanalytic(::Type{<:Pulse}) = IsAnalytic()
imanalytic(::Type{<:Pulse}) = IsAnalytic()
isprimitive(::Type{<:Pulse}) = IsPrimitive()

@inline intensity_point(p::Pulse, x,y) = κ(p::Pulse, x)*κ(p::Pulse, y)
@inline visibility_point(p::Pulse, u,v) = ω(p::Pulse, u)*ω(p::Pulse, u)

include("pulse.jl")
include("intensitymap.jl")
include("polarizedmap.jl")
