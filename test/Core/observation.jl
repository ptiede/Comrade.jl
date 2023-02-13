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
    dcoh = extract_coherency(obsavg)
    plot(dcoh)
    show(dcoh)

    @test getuv(arrayconfig(lcamp))[1] == getuv(arrayconfig(cphase))[1]
    @test getuv(arrayconfig(amp))[1] == getuv(arrayconfig(cphase))[1]
    @test getuv(arrayconfig(vis))[1] == getuv(arrayconfig(cphase))[1]


    m = stretched(Gaussian(), 1e-10,1e-10)
    plot(m, vis)
    plot(m, amp)
    plot(m, lcamp)
    plot(m,cphase)
    jt = jonesT(TransformCache(dcoh))
    mp = JonesModel(jt, PolarizedModel(m, ZeroModel(), ZeroModel(), 0.1*Gaussian()), TransformCache(dcoh))
    plot(mp, dcoh)
    #residual(m, vis)
    residual(m, amp)
    residual(m, cphase)
    residual(m, lcamp)
    residual(mp, dcoh)
    chi2(mp, dcoh)

    @test length(vis) == length(obsavg.data)
    @test length(amp) == length(obsavg.amp)
    @test length(extract_cphase(obsavg; count="min")) == length(obsavg.cphase)
    @test length(extract_lcamp(obsavg; count="min")) == length(obsavg.logcamp)
    #@test Set(stations(vis)) == Set(Symbol.(collect(get(obsavg.tarr, "site"))))
    @test mean(getdata(amp, :measurement)) == mean(get(obsavg.amp, :amp))

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
    @test visibilities(m, ac) ≈ visibilities(m, (U=u, V=v))

    @test visibility(m, ac.data[1]) ≈ visibility(m, (U=u[1], V=v[1]))


    u1 = cphase[:U1]
    v1 = cphase[:V1]
    u2 = cphase[:U2]
    v2 = cphase[:V2]
    u3 = cphase[:U3]
    v3 = cphase[:V3]
    bispectra(m, (U=u1, V=v1), (U=u2, V=v2), (U=u3, V=v3))

    @testset "RadioLikelihood" begin

        m = stretched(Gaussian(), μas2rad(20.0), μas2rad(20.0))
        f(_) = stretched(Gaussian(), μas2rad(20.0), μas2rad(20.0))
        lamp = RadioLikelihood(f, amp)
        lcphase = RadioLikelihood(f, cphase)
        llcamp = RadioLikelihood(f, lcamp)
        lclose1 = RadioLikelihood(f, cphase, lcamp)
        lclose2 = RadioLikelihood(f, lcamp, cphase)

        θ = NamedTuple{()}(())
        @test logdensityof(lclose1,θ) ≈ logdensityof(lclose2, θ)
        @test logdensityof(lclose1, θ ) ≈ logdensityof(lcphase, θ) + logdensityof(llcamp, θ)
        @inferred logdensityof(lclose1, θ)
        @inferred logdensityof(lclose2, θ)
        @inferred logdensityof(lamp, θ)

    end

end
