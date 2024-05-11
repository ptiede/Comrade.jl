using PythonCall
using Pyehtim


@testset "observation" begin

    println("ehtim:")

    obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "../test_data.uvfits"))

    obsavg = scan_average(obs)
    obsavg.add_cphase(snrcut=3, count="min")
    obsavg.add_amp()
    obsavg.add_logcamp(snrcut=3, count="min")

    obspol = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "../test_data.uvfits"); polrep="circ")
    obspolavg = scan_average(obspol)


    vis1 = Comrade.extract_vis(obsavg)
    vis  = extract_table(obsavg, ComplexVisibilities())
    @test vis[:measurement] ≈ vis1[:measurement]
    @test vis[:noise] ≈ vis1[:noise]
    @test vis[:baseline].U ≈ vis1[:baseline].U
    @test vis[:baseline].V ≈ vis1[:baseline].V
    plot(vis)
    show(vis)

    amp1 = Comrade.extract_amp(obsavg)
    amp = extract_table(obsavg, VisibilityAmplitudes())
    @test amp[:measurement] ≈ amp1[:measurement]
    @test amp[:noise] ≈ amp1[:noise]
    @test amp[:baseline].U ≈ amp1[:baseline].U
    @test amp[:baseline].V ≈ amp1[:baseline].V
    plot(amp)
    show(amp)

    cphase1 = Comrade.extract_cphase(obsavg; snrcut=3.0)
    cphase = extract_table(obsavg, ClosurePhases(;snrcut=3.0))
    @test cphase[:measurement] ≈ cphase1[:measurement]
    @test cphase[:noise] ≈ cphase1[:noise]
    @test cphase[:baseline].U1 ≈ cphase1[:baseline].U1
    @test cphase[:baseline].V1 ≈ cphase1[:baseline].V1
    @test cphase[:baseline].U2 ≈ cphase1[:baseline].U2
    @test cphase[:baseline].V2 ≈ cphase1[:baseline].V2
    @test cphase[:baseline].U3 ≈ cphase1[:baseline].U3
    @test cphase[:baseline].V3 ≈ cphase1[:baseline].V3
    plot(cphase)
    show(cphase)


    lcamp1 = Comrade.extract_lcamp(obsavg; snrcut=3.0)
    lcamp  = extract_table(obsavg, LogClosureAmplitudes(;snrcut=3.0))
    @test lcamp[:measurement] ≈ lcamp1[:measurement]
    @test lcamp[:noise] ≈ lcamp1[:noise]
    @test lcamp[:baseline].U1 ≈ lcamp1[:baseline].U1
    @test lcamp[:baseline].V1 ≈ lcamp1[:baseline].V1
    @test lcamp[:baseline].U2 ≈ lcamp1[:baseline].U2
    @test lcamp[:baseline].V2 ≈ lcamp1[:baseline].V2
    @test lcamp[:baseline].U3 ≈ lcamp1[:baseline].U3
    @test lcamp[:baseline].V3 ≈ lcamp1[:baseline].V3
    @test lcamp[:baseline].U4 ≈ lcamp1[:baseline].U4
    @test lcamp[:baseline].V4 ≈ lcamp1[:baseline].V4
    plot(lcamp)
    show(lcamp)

    dcoh1 = Comrade.extract_coherency(obspolavg)
    dcoh  = extract_table(obspolavg, Coherencies())
    @test dcoh[:measurement].:1 ≈ dcoh1[:measurement].:1
    @test dcoh[:noise].:1 ≈ dcoh1[:noise].:1
    @test dcoh[:baseline].U ≈ dcoh1[:baseline].U
    @test dcoh[:baseline].V ≈ dcoh1[:baseline].V

    plot(dcoh)
    show(dcoh)



    vis2, amp2, lcamp2, cphase2 = extract_table(obsavg, ComplexVisibilities(), VisibilityAmplitudes(),
                                                        ClosurePhases(; snrcut=3), LogClosureAmplitudes(;snrcut=3))


    @test getuv(arrayconfig(lcamp))[1] == getuv(arrayconfig(cphase))[1]
    @test getuv(arrayconfig(amp))[1] == getuv(arrayconfig(cphase))[1]
    @test getuv(arrayconfig(vis))[1] == getuv(arrayconfig(cphase))[1]


    m = stretched(Gaussian(), 1e-10,1e-10)
    plot(m, vis)
    plot(m, amp)
    plot(m, lcamp)
    plot(m,cphase)
    jt = jonesT(ResponseCache(dcoh))
    mp = Comrade.VLBIModel(JonesModel(jt), PolarizedModel(m, ZeroModel(), ZeroModel(), 0.1*Gaussian()))
    plot(mp, dcoh)
    #residual(m, vis)
    residual(m, amp)
    residual(m, cphase)
    residual(m, lcamp)
    residual(mp, dcoh)
    chi2(mp, dcoh)

    @test length(vis) == length(obsavg.data)
    @test length(amp) == length(obsavg.amp)
    # @test length(Comrade.extract_cphase(obsavg; snrcut=3, count="min")) == length(obsavg.cphase)
    # @test length(Comrade.extract_lcamp(obsavg; snrcut=3, count="min")) == length(obsavg.logcamp)
    #@test Set(sites(vis)) == Set(Symbol.(collect(get(obsavg.tarr, "site"))))
    # @test mean(getdata(amp, :measurement)) ≈ pyconvert(Float64, mean(obsavg.amp["amp"]))

    println("observation: ")
    uvpositions(vis[1])
    uvpositions(amp[1])
    uvpositions(cphase[1])
    uvpositions(lcamp[1])

    baselines(cphase[1])
    baselines(lcamp[1])

    #Test amplitude which is already debiased for very very dumb reasons
    @test sqrt(abs(visibility(vis[1]))^2 - vis[1].noise^2) ≈ amp[1].measurement
    ac = arrayconfig(vis)
    plot(ac)

    u,v = getuv(ac)
    @test visibilitymap(m, ac) ≈ visibilitymap(m, (U=u, V=v))

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
