using Distributions
using VLBIFiles
using VLBIImagePriors

_gauge_test_data() = load_test_table()

function _gauge_test_sky()
    tm(θ, meta) = θ.f1 * stretched(Gaussian(), θ.σ1, θ.σ1)
    prior = (f1 = VLBIUniform(0.5, 1.5), σ1 = VLBIUniform(μas2rad(1.0), μas2rad(100.0)))
    return SkyModel(tm, prior, spatialgrid(μas2rad(150.0), μas2rad(150.0), 32, 32))
end

# A station phase split into a track-long offset and a per-scan residual. Only the residual
# carries a reference, so the offsets share one level the data never see.
@instrument function gauge_split(; gauge = StationPhase(), offset_refant = NoReference())
    return @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)))
        gpμ ~ ArrayPrior(
            IIDSitePrior(TrackSeg(), DiagonalVonMises(0.0, inv(π^2)));
            refant = offset_refant, gauge
        )
        gp ~ ArrayPrior(
            IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2)); init = FixedInit(0.0));
            refant = SensitiveReference(0.0), gauge
        )
        return Jones(exp(complex(lg, gpμ + gp)))
    end
end

@testset "Phase gauge on an instrument model" begin
    dvis = _gauge_test_data()
    skym = _gauge_test_sky()
    t = dvis

    regauge(m, gaugefix) = InstrumentModel(m.jones; refbasis = m.refbasis, gaugefix)

    @testset "an unreferenced offset is reported" begin
        m = gauge_split()
        @test_throws "rank deficient by 1" bind_instrument(m, t)
        @test_throws "gpμ[site = AA" bind_instrument(m, t)
        @test_throws "gaugefix = PinGauge()" bind_instrument(m, t)
        # the analysis runs only on the terms that declare themselves
        @test bind_instrument(gauge_split(; gauge = nothing), t) isa Tuple
    end

    @testset "PinGauge identifies the phase" begin
        _, pinned = bind_instrument(regauge(gauge_split(), PinGauge()), t)
        _, loose = bind_instrument(gauge_split(; gauge = nothing), t)
        # a wrapped angle costs two flat coordinates, so one pin takes two off
        @test Comrade.dimension(asflat(pinned)) == Comrade.dimension(asflat(loose)) - 2

        # pinning the offset by hand instead reaches the same parameterization
        _, byhand = bind_instrument(
            gauge_split(; offset_refant = SensitiveReference(0.0)), t
        )
        @test Comrade.dimension(asflat(pinned)) == Comrade.dimension(asflat(byhand))

        # nothing is left over, and the pinned entry no longer varies
        names = (:gpμ, :gp)
        obs = NamedTuple{names}(map(n -> getproperty(pinned, n), names))
        @test isempty(Comrade.gauge_pins(Comrade.gauge_terms(obs, names), t))
        x, y = rand(pinned), rand(pinned)
        held = [i for i in eachindex(parent(x.gpμ)) if parent(x.gpμ)[i] == parent(y.gpμ)[i]]
        @test Comrade.sites(x.gpμ)[held] == [:AA]

        post = VLBIPosterior(skym, regauge(gauge_split(), PinGauge()), dvis)
        @test Comrade.dimension(asflat(post)) ==
            Comrade.dimension(asflat(VLBIPosterior(skym, gauge_split(; gauge = nothing), dvis))) - 2
    end

    @testset "gaugefix takes two values" begin
        m = gauge_split()
        @test_throws "must be ReportGauge() or PinGauge()" regauge(m, :pin)
        @test regauge(m, PinGauge()).gaugefix === PinGauge()
        @test m.gaugefix === ReportGauge()
    end

end
