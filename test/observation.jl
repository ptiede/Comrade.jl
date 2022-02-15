@testset "observation" begin

    println("ehtim:")

    load_ehtim()
    obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "test_data.uvfits"))
    obs.add_scans()
    obsavg = obs.avg_coherent(0.0, scan_avg=true)

    vis = extract_vis(obsavg)
    plot(vis)
    show(vis)
    amp = extract_amp(obsavg)
    plot(amp)
    show(amp)
    cphase = extract_cphase(obsavg)
    plot(cphase)
    show(cphase)
    lcamp = extract_lcamp(obsavg)
    plot(lcamp)
    show(lcamp)

    plot(stretched(Gaussian(), 1e-10, 1e-10), amp)
    residual(stretched(Gaussian(), 1e-10, 1e-10), amp)
    residual(stretched(Gaussian(), 1e-10, 1e-10), cphase)

    @test length(vis) == length(obsavg.data)
    @test length(amp) == length(obsavg.amp)
    @test length(cphase) == length(obsavg.cphase)
    @test length(lcamp) == length(obsavg.logcamp)
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




    @testset "RadioLikelihood" begin
        lamp = RadioLikelihood(amp)
        lcphase = RadioLikelihood(cphase)
        llcamp = RadioLikelihood(lcamp)
        lclose1 = RadioLikelihood(cphase, lcamp)
        lclose2 = RadioLikelihood(lcamp, cphase)

        m = stretched(Gaussian(), μas2rad(20.0), μas2rad(20.0))
        logdensity(lcphase, m)
        logdensity(llcamp, m)
        logdensity(lamp, m)

        @test logdensity(lclose1,m) ≈ logdensity(lclose2,m)
        @test logdensity(lclose1,m ) ≈ logdensity(lcphase, m) + logdensity(llcamp, m)

    end

end
