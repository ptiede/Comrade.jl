using Distributions
using Enzyme
using FiniteDifferences

# `t` as an observation of the source `name` at declination `dec`.
function with_source(t::VisTable, name; dec = info(t).source.dec)
    a, d = info(t), domain(t)
    c = Comrade.ComradeBase.coords(d)
    g = stations(t)
    b = ArrayInfo(; a.stations, a.channels, a.scans, a.tracks, source = Source(name, a.source.ra, dec), a.epoch)
    return VisTable(
        t.pol, Comrade._dense(measurement(t)); weight = Comrade._dense(weight(t)),
        u = c.u, v = c.v, Ti = d.Ti, ant1 = d.ant1, ant2 = d.ant2, info = b,
        elevation = parent(g.elevation), parallactic = parent(g.parallactic),
    )
end

@testset "Fields" begin
    _, dvis, = load_data()
    skym = SkyModel(test_model, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 64, 64))
    gainterm(seg) = Jones(
        (x, c) -> exp(complex(x.lg, x.gp));
        lg = ArrayPrior(IIDSitePrior(seg, Normal(0.0, 0.1))),
        gp = ArrayPrior(IIDSitePrior(seg, DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
    )
    cal = with_source(dvis, :CAL)

    @testset "simultaneous fields share the instrument" begin
        intm = InstrumentModel(gainterm(ScanSeg()))
        post1 = VLBIPosterior(skym, intm, dvis)
        post2 = VLBIPosterior(intm, Field(:M87, skym, dvis), Field(:CAL, skym, cal))
        θ1 = prior_sample(post1)
        θ2 = (; sky = (; M87 = θ1.sky, CAL = θ1.sky), θ1.instrument)
        @test dimension(post2) == dimension(post1) + length(θ1.sky)
        @test loglikelihood(post2, θ2) ≈ 2 * loglikelihood(post1, θ1)
        @test logdensityof(post2, θ2) ≈ 2 * loglikelihood(post1, θ1) + Comrade.logprior(post2, θ2)
        @test keys(prior_sample(post2).sky) == (:M87, :CAL)

        @test keys(simulate_observation(post2, θ2; add_thermal_noise = false)) == (:M87, :CAL)
        r = Comrade.residuals(post2, θ2)
        @test parent(measurement(only(r.CAL))) ≈ parent(measurement(only(Comrade.residuals(post1, θ1))))
        @test collect(chi2(post2, θ2).M87) ≈ collect(chi2(post1, θ1))
        @test keys(instrumentmodel(post2, θ2)) == (:M87, :CAL)

        tpost = asflat(post2)
        x = prior_sample(tpost)
        gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(tpost), x)
        gfd, = FiniteDifferences.grad(central_fdm(5, 1), tpost, x)
        @test gz ≈ gfd
    end

    @testset "fields at different times" begin
        late = later(cal, 1)
        for (seg, n) in ((ScanSeg(), 2), (EpochSeg(), 2), (TrackSeg(), 1))
            intm = InstrumentModel(gainterm(seg))
            single = dimension(VLBIPosterior(skym, intm, dvis)) - length(test_prior())
            post = VLBIPosterior(intm, Field(:M87, skym, dvis), Field(:CAL, skym, late))
            @test dimension(post) - 2 * length(test_prior()) == n * single
        end
        post = VLBIPosterior(InstrumentModel(gainterm(ScanSeg())), Field(:M87, skym, dvis), Field(:CAL, skym, late))
        L = post.instrumentmodel.fields
        @test L.CAL.layout.coords.Ti ≈ L.M87.layout.coords.Ti .+ 24
    end

    @testset "direction-dependent terms" begin
        dec0 = deg2rad(info(dvis).source.dec)
        A = Jones((x, c) -> exp(x.a * (c.dec - dec0)); a = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 1.0))))
        north = with_source(dvis, :N; dec = info(dvis).source.dec + 1)
        post = VLBIPosterior(InstrumentModel(A), Field(:M87, skym, dvis), Field(:N, skym, north))
        J = post.instrumentmodel.fields
        @test all(≈(dec0 + deg2rad(1)), J.N.layout.coords.dec)
        x = prior_sample(post).instrument
        parent(x.a) .= 0.5
        vis = parent(measurement(dvis))
        @test Comrade.apply_instrument(vis, J.M87, (; instrument = x)) ≈ vis
        @test Comrade.apply_instrument(vis, J.N, (; instrument = x)) ≈ vis .* exp(2 * 0.5 * deg2rad(1))
    end

    @testset "errors" begin
        intm = InstrumentModel(gainterm(ScanSeg()))
        @test_throws "distinct names" VLBIPosterior(intm, Field(:A, skym, dvis), Field(:A, skym, cal))
        @test_throws "no data products" Field(:A, skym)
        @test_throws "one table of visibilities" VLBIPosterior(intm, Field(:A, skym, dvis, amplitudes(later(dvis, 1))))
        @test_throws "several fields" VLBIPosterior(InstrumentModel(gainterm(ScanSeg()); selfcal = true), Field(:A, skym, dvis), Field(:B, skym, cal))
    end
end
