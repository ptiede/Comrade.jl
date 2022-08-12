@testset "observation" begin

    println("ehtim:")

    load_ehtim()
    obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "../test_data.uvfits"))
    obs.add_scans()
    obsavg = obs.avg_coherent(0.0, scan_avg=true)
    obsavg.add_cphase(count="min")
    obsavg.add_amp()
    obsavg.add_logcamp(count="min")

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

    @test getuv(arrayconfig(lcamp))[1] == getuv(arrayconfig(cphase))[1]
    @test getuv(arrayconfig(amp))[1] == getuv(arrayconfig(cphase))[1]
    @test getuv(arrayconfig(vis))[1] == getuv(arrayconfig(cphase))[1]


    m = stretched(Gaussian(), 1e-10,1e-10)
    plot(m, vis)
    plot(m, amp)
    plot(m, lcamp)
    plot(m,cphase)
    #residual(m, vis)
    residual(m, amp)
    residual(m, cphase)
    residual(m, lcamp)

    @test length(vis) == length(obsavg.data)
    @test length(amp) == length(obsavg.amp)
    @test length(extract_cphase(obsavg; count="min")) == length(obsavg.cphase)
    @test length(extract_lcamp(obsavg; count="min")) == length(obsavg.logcamp)
    #@test Set(stations(vis)) == Set(Symbol.(collect(get(obsavg.tarr, "site"))))
    @test mean(getdata(amp, :amp)) == mean(get(obsavg.amp, :amp))

    println("observation: ")
    uvpositions(vis[1])
    uvpositions(amp[1])
    uvpositions(cphase[1])
    uvpositions(lcamp[1])

    baselines(cphase[1])
    baselines(lcamp[1])

    #Test amplitude which is already debiased for very very dumb reasons
    @test sqrt(abs(visibility(vis[1]))^2 - vis[1].error^2) ≈ amplitude(amp[1])
    ac = arrayconfig(vis)
    plot(ac)

    u,v = getuv(ac)
    @test visibilities(m, ac) ≈ visibilities(m, u, v)

    @test visibility(m, ac.data[1]) ≈ visibility(m, u[1], v[1])


    u1 = cphase[:u1]
    v1 = cphase[:v1]
    u2 = cphase[:u2]
    v2 = cphase[:v2]
    u3 = cphase[:u3]
    v3 = cphase[:v3]
    bispectra(m, u1, v1, u2, v2, u3, v3)

    @testset "RadioLikelihood" begin

        m = stretched(Gaussian(), μas2rad(20.0), μas2rad(20.0))
        logdensityof(lcphase, m)
        logdensityof(llcamp, m)
        logdensityof(lamp, m)

        lamp = RadioLikelihood(m, amp)
        lcphase = RadioLikelihood(m, cphase)
        llcamp = RadioLikelihood(m, lcamp)
        lclose1 = RadioLikelihood(m, cphase, lcamp)
        lclose2 = RadioLikelihood(m, lcamp, cphase)


        @test logdensityof(lclose1,m) ≈ logdensityof(lclose2,m)
        @test logdensityof(lclose1,m ) ≈ logdensityof(lcphase, m) + logdensityof(llcamp, m)
        @inferred logdensityof(lclose1, m)
        @inferred logdensityof(lclose2, m)
        @inferred logdensityof(lamp, m)

    end

end
