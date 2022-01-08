@testset "observation" begin

    println("ehtim:")

    load_ehtim()
    obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "test_data.uvfits"))
    obs.add_scans()
    obsavg = obs.avg_coherent(0.0, scan_avg=true)

    vis = extract_vis(obsavg)
    show(vis)
    amp = extract_amp(obsavg)
    show(amp)
    cphase = extract_cphase(obsavg)
    show(cphase)
    lcamp = extract_lcamp(obsavg)
    show(lcamp)

    @test nsamples(vis) == length(obsavg.data)
    @test nsamples(amp) == length(obsavg.amp)
    @test nsamples(cphase) == length(obsavg.cphase)
    @test nsamples(lcamp) == length(obsavg.logcamp)
    @test Set(stations(vis)) == Set(Symbol.(collect(get(obsavg.tarr, "site"))))
    @test mean(getdata(amp, :amp)) == mean(get(obsavg.amp, :amp))

    println("observation: ")
    uvpositions(vis[1])
    uvpositions(amp[1])
    uvpositions(cphase[1])
    uvpositions(lcamp[1])

    baselines(cphase[1])
    baselines(lcamp[1])

    #Test amplitude which is already debiased for very very dumb reasons
    @test sqrt(abs(visibility(vis[1]))^2 - vis[1].error^2) == amplitude(amp[1])
    ac = arrayconfig(vis)
end
