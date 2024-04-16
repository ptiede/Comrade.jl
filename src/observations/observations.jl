
using AstroTime: modified_julian

export uvpositions, sites, getdata, arrayconfig,
       getuv, baselines, scantable, beamsize


include("dataproducts.jl")
include("datums.jl")
include("abstract.jl")
include("array.jl")
include("obstable.jl")
include("scantable.jl")











# # internal method that checks whether the triangle is closes
# function checktriangle(D1::EHTComplexVisibilityDatum,
#                        D2::EHTComplexVisibilityDatum,
#                        D3::EHTComplexVisibilityDatum)
#     b1 = D1.baseline
#     b2 = D2.baseline
#     b3 = D3.baseline
#     l = length(unique([b1..., b2..., b3...]))
#     @assert l == 3 "For a valid closure phase you need 3 unique sites not $l"
#     @assert (D1.time == D2.time == D3.time) "For a valid closure phase the times need to match"

# end


# """
#     bispectrum(d1::T, d2::T, d3::T) where {T<:EHTComplexVisibilityDatum}

# Finds the bispectrum of three visibilities. We will assume these form closed triangles,
# i.e. the phase of the bispectrum is a closure phase.
# """
# @inline function bispectrum(D1::EHTComplexVisibilityDatum, D2::EHTComplexVisibilityDatum, D3::EHTComplexVisibilityDatum)
#     checktriangle(D1, D2, D3)
#     visibility(D1)*visibility(D2)*visibility(D3)
# end






# """
#     closure_phase(D1::EHTComplexVisibilityDatum,
#                   D2::EHTComplexVisibilityDatum,
#                   D3::EHTComplexVisibilityDatum
#                   )

# Computes the closure phase of the three visibility datums.

# # Notes
# We currently use the high SNR Gaussian error approximation for the closure phase.
# In the future we may use the moment matching from Monte Carlo sampling.
# """
# function closure_phase(D1::EHTComplexVisibilityDatum,
#                       D2::EHTComplexVisibilityDatum,
#                       D3::EHTComplexVisibilityDatum)

#     checktriangle(D1,D2,D3)

#     amp1 = amplitude(D1).amp
#     amp2 = amplitude(D2).amp
#     amp3 = amplitude(D3).amp
#     u1,v1 = uvpositions(D1)
#     u2,v2 = uvpositions(D2)
#     u3,v3 = uvpositions(D3)

#     bis = bispectrum(D1, D2, D3)
#     s12 = unique([D1.baseline..., D2.baseline...])
#     s123 = unique([s12..., D3.baseline...])
#     #Use the Gaussian approximation TODO hook this into Measurements.jl?
#     err = sqrt((D1.error/amp1)^2 + (D2.error/amp2)^2 + (D3.error/amp3)^2)
#     return EHTClosurePhaseDatum(angle(bis), err,
#                                 u1, v1, u2, v2, u3, v3,
#                                 time, s123)
# end







function _arrayconfig(data, angles, tarr, scans, bandwidth, ra, dec, mjd, source)
    u = getproperty(data, :U)
    v = getproperty(data, :V)
    times = getproperty(data, :T)
    error = getproperty(data, :error)
    baseline = getproperty(data, :baseline)
    frequency = getproperty(data, :F)
    uvsamples = StructArray{EHTArrayBaselineDatum}(T=times,
                                        U=u,
                                        V=v,
                                        F = frequency,
                                        baseline=baseline,
                                        error=error,
                                        elevation = StructArray(angles[1]),
                                        parallactic  = StructArray(angles[2])
                                    )
    return EHTArrayConfiguration(bandwidth, tarr, scans, mjd, ra, dec, source, :UTC, uvsamples)
end
