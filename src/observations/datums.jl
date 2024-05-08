export measurement, noise, baseline

"""
    AbstractVisibilityDatum
An abstract type for all VLBI data types. See [`Comrade.EHTVisibilityDatum`](@ref) for an example.
"""
abstract type AbstractVisibilityDatum{T} end
baseline(p::AbstractVisibilityDatum) = getfield(p, :baseline)
measurement(p::AbstractVisibilityDatum) = getfield(p, :measurement)
noise(p::AbstractVisibilityDatum) = getfield(p, :noise)
VLBISkyModels.polarization(p::AbstractVisibilityDatum) = getfield(p, :polarization)

# function Base.propertynames(p::AbstractVisibilityDatum)
#     return (propertynames(baseline(p))..., :measurement, :noise)
# end
# function Base.getproperty(p::AbstractVisibilityDatum, s::Symbol)
#     s == :measurement && return measurement(p)
#     s == :noise       && return noise(p)
#     return getproperty(baseline(p), s)
# end

build_datum(F::Type{<:AbstractVisibilityDatum}, m, e, b) = F(m, e, b)

abstract type AbstractSinglePolDatum{P,S} <: AbstractVisibilityDatum{S} end
abstract type ClosureProducts{P,T} <: AbstractSinglePolDatum{P,T} end

abstract type AbstractBaselineDatum end

"""
    $(TYPEDEF)

A Datum for a single coherency matrix

# Fields
$(FIELDS)

"""
Base.@kwdef struct EHTCoherencyDatum{S, B<:AbstractBaselineDatum, M<:SMatrix{2,2,Complex{S}}, E<:SMatrix{2,2,S}} <: Comrade.AbstractVisibilityDatum{S}
    """
    coherency matrix, with entries in Jy
    """
    measurement::M
    """
    visibility uncertainty matrix, with entries in Jy
    """
    noise::E
    """
    baseline information
    """
    baseline::B
end
VLBISkyModels.polarization(b::EHTCoherencyDatum) = b.baseline.polbasis


"""
    $(TYPEDEF)

A struct holding the information for a single measured complex visibility.

# FIELDS
$(FIELDS)

"""
Base.@kwdef struct EHTVisibilityDatum{Pol, S<:Number, B<:AbstractBaselineDatum} <: AbstractSinglePolDatum{Pol,S}
    """
    Complex Vis. measurement (Jy)
    """
    measurement::Complex{S}
    """
    noise of the complex vis (Jy)
    """
    noise::S
    """
    baseline information
    """
    baseline::B
end
EHTVisibilityDatum(m::Complex{S}, n, b) where {S} = EHTVisibilityDatum{:I, S, typeof(b)}(m, n, b)




"""
    $(TYPEDEF)

A struct holding the information for a single measured visibility amplitude.

# FIELDS
$(FIELDS)

"""
Base.@kwdef struct EHTVisibilityAmplitudeDatum{P, S<:Number, B<:AbstractBaselineDatum} <: AbstractSinglePolDatum{P,S}
    """
    amplitude (Jy)
    """
    measurement::S
    """
    noise of the visibility amplitude (Jy)
    """
    noise::S
    """
    baseline information
    """
    baseline::B
end
EHTVisibilityAmplitudeDatum(m::Number, n::Number, b) = EHTVisibilityAmplitudeDatum{:I, typeof(m), typeof(b)}(m, n, b)


"""
    $(TYPEDEF)

A Datum for a single closure phase.

# Fields
$(FIELDS)

"""
Base.@kwdef struct EHTClosurePhaseDatum{P,S<:Number, B<:AbstractBaselineDatum} <: ClosureProducts{P,S}
    """
    closure phase (rad)
    """
    measurement::S
    """
    noise of the closure phase assuming the high-snr limit
    """
    noise::S
    """
    baselines for the closure phase
    """
    baseline::NTuple{3, B}
end
EHTClosurePhaseDatum(m::Number, n::Number, b) = EHTClosurePhaseDatum{:I, typeof(m), typeof(b)}(m, n, b)


"""
    triangle(b::EHTClosurePhaseDatum)

Returns the sites used in the closure phase triangle.
"""
triangle(b::EHTClosurePhaseDatum) = map(x->first(getproperty(x, :baseline)), baseline(b))



"""
    $(TYPEDEF)

A Datum for a single log closure amplitude.

# $(FIELDS)

"""
Base.@kwdef struct EHTLogClosureAmplitudeDatum{P, S<:Number, B<:AbstractBaselineDatum} <: ClosureProducts{P,S}
    """
    log-closure amplitude
    """
    measurement::S
    """
    log-closure amplitude noise in the high-snr limit
    """
    noise::S
    """
    baselines for the closure phase
    """
    baseline::NTuple{4, B}
end
EHTLogClosureAmplitudeDatum(m::Number, n::Number, b) = EHTLogClosureAmplitudeDatum{:I, typeof(m), typeof(b)}(m, n, b)
