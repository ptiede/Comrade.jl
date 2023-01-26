"""
    $(TYPEDEF)
An abstract type for all VLBI interfermetry data types. See [EHTVisibilityDatum](@ref) for an example.
"""
abstract type AbstractInterferometryDatum{T} end

measurement(d::AbstractInterferometryDatum) = d.measurement
error(d::AbstractInterferometryDatum) = d.error
baselinedatum(d::AbstractInterferometryDatum) = d.baselinedatum

abstract type AbstractVisibilityDatum{T} <: AbstractInterferometryDatum{T} end
abstract type AbstractLinearPolDatum{S<:AbstractVisibilityDatum, T} <: AbstractInterferometryDatum{T} end
abstract type AbstractCrossPolDatum{S,T} <: AbstractInterferometryDatum{T} end

abstract type ClosureProducts{T} <: AbstractInterferometryDatum{T} end

abstract type BaselineDatum end

"""
    $(TYPEDEF)

A single datum of an `ArrayConfiguration`
"""
struct EHTBaselineDatum{T,P,E,V} <: BaselineDatum
    """
    u position of the data point in λ
    """
    U::T
    """
    v position of the data point in λ
    """
    V::T
    """
    time of the data point in (Hr)
    """
    T::T
    """
    frequency of the data point (Hz)
    """
    F::T
    """
    Station codes of the baseline (u,v)
    """
    baseline::Tuple{Symbol, Symbol}
    """
    Polarization basis of the measurement
    """
    polbasis::P
    """
    elevation of baselines
    """
    elevation::Tuple{V,V}
    """
    parallactic angle of baslines
    """
    parallactic::Tuple{V,V}
    function EHTBaselineDatum(u, v, time, freq, baseline, polbasis, elevation, parallactic)
        tt, ft, ut, vt = promote(time, freq, u, v)
        T = typeof(tt)
        V = typeof(elevation[1])
        E = typeof(error)
        return new{T,E,V}(ut, vt, tt, ft, baseline, polbasis, elevation, parallactic)
    end
end

uvpositions(b::BaselineDatum) = (U=b.U, V=b.V)



"""
    $(TYPEDEF)

A struct holding the information for a single measured visibility.

# $(FIELDS)

"""
Base.@kwdef struct EHTVisibilityDatum{S<:Number,D} <: AbstractVisibilityDatum{S}
    """
    real component of the visibility (Jy)
    """
    measurement::Complex{S}
    """
    error of the visibility (Jy)
    """
    error::S
    """
    Baseline datum
    """
    baselinedatum::P
end

function build(S::Type{<:AbstractInterferometryDatum}, b::ArrayBaselineDatum, measurement, error)
    S(measurement = measurement,
      error = error,
      baselinedatum = b)
end


"""
    uvpositions(datum::AbstractVisibilityDatum)

Get the uvp positions of an inferometric datum.
"""
uvpositions(D::AbstractVisibilityDatum) = uvpositions(baselinedatum(D))


"""
    $(TYPEDEF)

A struct holding the information for a single measured visibility amplitude.

# FIELDS
$(FIELDS)

"""
Base.@kwdef struct EHTVisibilityAmplitudeDatum{S<:Number,D} <: AbstractVisibilityDatum{S}
    """
    amplitude (Jy)
    """
    measurement::S
    """
    error of the visibility amplitude (Jy)
    """
    error::S
    """
    baseline datum
    """
    baselinedatum::D
end


"""
    $(TYPEDEF)

A Datum for a single coherency matrix

# Fields
$(FIELDS)

"""
Base.@kwdef struct EHTCoherencyDatum{S,M<:SMatrix{2,2,Complex{S}}, E<:SMatrix{2,2,S},D} <: Comrade.AbstractInterferometryDatum{S}
    """
    coherency matrix, with entries in Jy
    """
    measurement::M
    """
    visibility uncertainty matrix, with entries in Jy
    """
    error::E
    """
    baseline datum
    """
    baselinedatum::D
end

"""
    $(TYPEDEF)

A Datum for a single log closure amplitude.

# $(FIELDS)

"""
Base.@kwdef struct EHTLogClosureAmplitudeDatum{S<:Number,P,D} <: ClosureProducts{S}
    """
    log-closure amplitude
    """
    measurement::S
    """
    log-closure amplitude error in the high-snr limit
    """
    error::S
    """
    baseline datum for the first station
    """
    baseline1::D
    """
    baseline datum for the second station
    """
    baseline2::D
    """
    baseline datum for the third station
    """
    baseline3::D
    """
    baseline datum for the fourth station
    """
    baseline4::D
end

"""
    baselines(CP::EHTLogClosureAmplitudeDatum)

Returns the baselines used for a single closure phase datum
"""
function baselines(CP::EHTLogClosureAmplitudeDatum)
    return (baselines(CP.baseline1),
            baselines(CP.baseline2),
            baselines(CP.baseline3),
            baselines(CP.baseline4))
end

function uvpositions(datum::EHTLogClosureAmplitudeDatum)
    (uvpositions(datum.baseline1),
     uvpositions(datum.baseline2),
     uvpositions(datum.baseline3),
     uvpositions(datum.baseline4))
end


"""
    $(TYPEDEF)

A Datum for a single closure phase.

# Fields
$(FIELDS)

"""
Base.@kwdef struct EHTClosurePhaseDatum{S<:Number, D} <: ClosureProducts{S}
    """
    closure phase (rad)
    """
    measurement::S
    """
    error of the closure phase assuming the high-snr limit
    """
    error::S
    """
    baseline datum for the first station
    """
    baseline1::D
    """
    baseline datum for the second station
    """
    baseline2::D
    """
    baseline datum for the third station
    """
    baseline3::D
end

# internal method that checks whether the triangle is closes
function checktriangle(
    D1::EHTVisibilityDatum,
    D2::EHTVisibilityDatum,
    D3::EHTVisibilityDatum)
    b1 = D1.baseline
    b2 = D2.baseline
    b3 = D3.baseline
    l = length(unique([b1..., b2..., b3...]))
    @assert l == 3 "For a valid closure phase you need 3 unique stations not $l"
    @assert (D1.time == D2.time == D3.time) "For a valid closure phase the times need to match"
    return nothing
end



"""
    baselines(CP::EHTClosurePhaseDatum)

Returns the baselines used for a single closure phase datum
"""
function baselines(CP::EHTClosurePhaseDatum)
    return (baselines(CP.baseline1),
            baselines(CP.baseline2),
            baselines(CP.baseline3))
end

function uvpositions(datum::EHTClosurePhaseDatum)
    return (uvpositions(datum.baseline1),
            uvpositions(datum.baseline2),
            uvpositions(datum.baseline3),)
end

"""
    bispectrum(d1::T, d2::T, d3::T) where {T<:EHTVisibilityDatum}

Finds the bispectrum of three visibilities.
"""
@inline function bispectrum(D1::EHTVisibilityDatum, D2::EHTVisibilityDatum, D3::EHTVisibilityDatum)
    checktriangle(D1, D2, D3)
    visibility(D1)*visibility(D2)*visibility(D3)
end

"""
    closure_phase(D1::EHTVisibilityDatum,
                  D2::EHTVisibilityDatum,
                  D3::EHTVisibilityDatum
                  )

Computes the closure phase of the three visibility datums.

# Notes
We currently use the high SNR Gaussian error approximation for the closure phase.
In the future we may use the moment matching from Monte Carlo sampling.
"""
function closure_phase(D1::EHTVisibilityDatum,
                      D2::EHTVisibilityDatum,
                      D3::EHTVisibilityDatum)

    checktriangle(D1,D2,D3)

    amp1 = measurement(amplitude(D1))
    amp2 = measurement(amplitude(D2))
    amp3 = measurement(amplitude(D3))

    bis = bispectrum(D1, D2, D3)
    #Use the Gaussian approximation TODO hook this into Measurements.jl?
    err = sqrt((D1.error/amp1)^2 + (D2.error/amp2)^2 + (D3.error/amp3)^2)
    return EHTClosurePhaseDatum(angle(bis), err,
                                baselinedatum(D1),
                                baselinedatum(D2),
                                baselinedatum(D3))
end
