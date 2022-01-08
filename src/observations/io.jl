


"""
$(SIGNATURES)
Load a ThemisPy style ascii EHT observation file.
"""
function load_tpy(file)
    data = readdlm(file, skipstart=1)
    bs1 = Symbol.(getindex.(data[:,5], Ref(1:2)))
    bs2 = Symbol.(getindex.(data[:,5], Ref(3:4)))
    baselines = tuple.(bs1, bs2)
    edata = StructArray{EHTVisibilityDatum{Float64}}(
                visr=float.(data[:,8]),
                visi=float.(data[:,10]),
                error=float.(data[:,9]),
                u=float.(data[:,6])*1e6,
                v=float.(data[:,7])*1e6,
                time=float.(data[:,4]),
                frequency=fill(227e9, size(data,1)),
                bandwidth=fill(4e6, size(data,1)),
                baseline=baselines
            )
    mjd =  Int(modified_julian(UTCEpoch(Int(data[1,2]), 1,1,0)).Δt+float(data[1,3]))
    return EHTObservation(data=edata, mjd=mjd, ra=180.0, dec=0.0, source=Symbol(data[1,1]))
end
