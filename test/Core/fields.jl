using Distributions
using Enzyme
using FiniteDifferences

# `t` as an observation of the source `name` at declination `dec`.
function with_source(t::VisTable, name; dec = info(t).source.dec)
    a = info(t)
    b = ArrayInfo(; a.stations, a.channels, a.scans, a.tracks, source = Source(name, a.source.ra, dec), a.epoch)
    return VisTable(t.pol, records(t), collect(baseimage(measurement(t))), collect(baseimage(weight(t))); info = b)
end

@testset "Fields" begin
    _, dvis, = load_data()
    skym = SkyModel(test_model, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 64, 64))
    gainterm(seg) = Jones(
        (x, c) -> exp(complex(x.lg, x.gp));
        lg = ArrayPrior(IIDStationPrior(seg, Normal(0.0, 0.1))),
        gp = ArrayPrior(IIDStationPrior(seg, DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
    )
    cal = with_source(dvis, :CAL)

    @testset "simultaneous fields share the instrument" begin
        intm = InstrumentModel(gainterm(ScanSeg()))
        post1 = VLBIPosterior(dvis; sky = skym, instrument = intm)
        post2 = VLBIPosterior((; M87 = dvis, CAL = cal); sky = (; M87 = skym, CAL = skym), instrument = intm)
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
            single = dimension(VLBIPosterior(dvis; sky = skym, instrument = intm)) - length(test_prior())
            post = VLBIPosterior((; M87 = dvis, CAL = late); sky = (; M87 = skym, CAL = skym), instrument = intm)
            @test dimension(post) - 2 * length(test_prior()) == n * single
        end
        post = VLBIPosterior((; M87 = dvis, CAL = late); sky = (; M87 = skym, CAL = skym), instrument = InstrumentModel(gainterm(ScanSeg())))
        L = post.instrumentmodel.fields
        @test L.CAL.layout.coords.Ti ≈ L.M87.layout.coords.Ti .+ 24
    end

    @testset "direction-dependent terms" begin
        dec0 = deg2rad(info(dvis).source.dec)
        A = Jones((x, c) -> exp(x.a * (c.dec - dec0)); a = ArrayPrior(IIDStationPrior(TrackSeg(), Normal(0.0, 1.0))))
        north = with_source(dvis, :N; dec = info(dvis).source.dec + 1)
        post = VLBIPosterior((; M87 = dvis, N = north); sky = (; M87 = skym, N = skym), instrument = InstrumentModel(A))
        J = post.instrumentmodel.fields
        @test all(≈(dec0 + deg2rad(1)), J.N.layout.coords.dec)
        x = prior_sample(post).instrument
        parent(x.a) .= 0.5
        vis = parent(measurement(dvis))
        @test Comrade.apply_instrument(vis, J.M87, x) ≈ vis
        @test Comrade.apply_instrument(vis, J.N, x) ≈ vis .* exp(2 * 0.5 * deg2rad(1))
    end

    @testset "several tables of one sky" begin
        intm = InstrumentModel(gainterm(ScanSeg()))
        late = later(dvis, 1)
        single = VLBIPosterior(dvis; sky = skym, instrument = intm)
        post = VLBIPosterior((; A = dvis, B = late); sky = skym, instrument = intm)
        θ = prior_sample(post)
        @test keys(θ.sky) == keys(prior_sample(single).sky)
        @test dimension(post) == length(test_prior()) + 2 * (dimension(single) - length(test_prior()))
        ideal = VLBIPosterior((; A = dvis, B = late); sky = skym)
        θi = prior_sample(ideal)
        @test loglikelihood(ideal, θi) ≈
            loglikelihood(VLBIPosterior(dvis; sky = skym), θi) + loglikelihood(VLBIPosterior(late; sky = skym), θi)
        @test keys(instrumentmodel(post, θ)) == (:A, :B)
        @test keys(simulate_observation(post, θ; add_thermal_noise = false)) == (:A, :B)
        @test keys(chi2(post, θ)) == (:A, :B)
        @test typeof(skymodel(post, θ)) == typeof(skymodel(single, θ))

        tpost = asflat(post)
        x = prior_sample(tpost)
        gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(tpost), x)
        gfd, = FiniteDifferences.grad(central_fdm(5, 1), tpost, x)
        @test gz ≈ gfd
        @test_throws "needs at least one table" VLBIPosterior((;); sky = skym)

        # Tables whose receivers differ at one station bind one instrument, though they cannot be joined.
        a = info(late)
        s1 = first(a.stations)
        r1 = only(receivers(s1))
        moved = Station(s1.name, s1.position, Receiver(MountEquatorial(), r1.feeds; r1.orientations))
        b = ArrayInfo(; stations = [moved; a.stations[2:end]], a.channels, a.scans, a.tracks, a.source, a.epoch)
        remounted = VisTable(late.pol, records(late), collect(baseimage(measurement(late))), collect(baseimage(weight(late))); info = b)
        two = VLBIPosterior((; A = dvis, B = remounted); sky = skym, instrument = intm)
        @test dimension(two) == dimension(post)
        @test logdensityof(two, θ) ≈ logdensityof(post, θ)
        @test_throws "overlapping bands" vcat(dvis, remounted)
    end

    @testset "frequency phase transfer across tables" begin
        ν0 = only(info(dvis).channels.frequency)
        dense = with_channel(dvis, ν0 / 2)
        lo = dense[Fr = At(ν0 / 2)]
        P = Jones(
            (x, c) -> exp(1im * (x.φ * c.Fr / ν0 + x.δ * ν0 / c.Fr + x.ψ));
            φ = ArrayPrior(IIDStationPrior(ScanSeg(), Normal(0.0, 1.0)); refant = SensitiveReference(0.0)),
            δ = ArrayPrior(IIDStationPrior(ScanSeg(), Normal(0.0, 0.1)); refant = SensitiveReference(0.0)),
            ψ = ArrayPrior(IIDStationPrior(TrackSeg() ⊗ BandSeg(bands(dense)...), Normal(0.0, 0.1)); refant = SensitiveReference(0.0)),
        )
        intm = InstrumentModel(P)
        post = VLBIPosterior((; lo, hi = dvis); sky = skym, instrument = intm)
        joined = VLBIPosterior(dense; sky = skym, instrument = intm)
        @test dimension(post) == dimension(joined)
        θ = prior_sample(post)
        @test logdensityof(post, θ) ≈ logdensityof(joined, θ)

        x = θ.instrument
        foreach(p -> parent(p) .= 0, (x.φ, x.δ, x.ψ))
        s = first(stations(dvis))
        xφ, xδ, xψ = x.φ[Ant = At(s)], x.δ[Ant = At(s)], x.ψ[Ant = At(s), Fr = At(ν0 / 2)]
        xφ .= 0.3
        xδ .= 0.05
        xψ .= -0.2
        for (t, J, θs) in ((lo, post.instrumentmodel.fields.lo, 0.3 / 2 + 0.05 * 2 - 0.2), (dvis, post.instrumentmodel.fields.hi, 0.3 + 0.05))
            vis = parent(measurement(t))
            side = map(b -> (b[1] == s) - (b[2] == s), baselines(t))
            @test Comrade.apply_instrument(vis, J, x) ≈ vis .* exp.(1im .* θs .* side)
        end
    end

    @testset "errors" begin
        intm = InstrumentModel(gainterm(ScanSeg()))
        @test_throws "visibility-domain data product" VLBIPosterior((; A = ()); sky = (; A = skym), instrument = intm)
        @test_throws "one table of visibilities" VLBIPosterior((; A = (dvis, amplitudes(later(dvis, 1)))); sky = (; A = skym), instrument = intm)
    end
end
