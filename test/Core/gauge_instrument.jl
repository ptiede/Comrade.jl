using Distributions
using VLBIFiles
using VLBIImagePriors

function _gauge_test_data()
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    return extract_table(uvd, Visibilities(; time_average = VLBI.GapBasedScans()))
end

function _gauge_test_sky()
    tm(θ, meta) = θ.f1 * stretched(Gaussian(), θ.σ1, θ.σ1)
    prior = (f1 = VLBIUniform(0.5, 1.5), σ1 = VLBIUniform(μas2rad(1.0), μas2rad(100.0)))
    return SkyModel(tm, prior, imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32))
end

# A station phase split into a track-long offset and a per-scan residual. Only the residual
# carries a reference, so the offsets share one level the data never see.
@instrument function gauge_split(; gauge = :phase, offset_refant = NoReference())
    return @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)))
        gpμ ~ ArrayPrior(
            IIDSitePrior(TrackSeg(), DiagonalVonMises(0.0, inv(π^2)));
            refant = offset_refant, gauge
        )
        gp ~ ArrayPrior(
            IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2)); init = FixedInit(0.0));
            refant = SEFDReference(0.0), gauge
        )
        return SingleStokesGain(exp(complex(lg, gpμ + gp)))
    end
end

@testset "Phase gauge on an instrument model" begin
    dvis = _gauge_test_data()
    skym = _gauge_test_sky()
    arr = arrayconfig(dvis)

    regauge(m, gaugefix) = InstrumentModel(m.jones, m.prior; refbasis = m.refbasis, gaugefix)

    @testset "an unreferenced offset is reported" begin
        m = gauge_split()
        @test_throws "rank deficient by 1" Comrade.set_array(m, arr)
        @test_throws "gpμ[site = AA" Comrade.set_array(m, arr)
        @test_throws "gaugefix = :pin" Comrade.set_array(m, arr)
        # the analysis runs only on the terms that declare themselves
        @test Comrade.set_array(gauge_split(; gauge = :none), arr) isa Tuple
    end

    @testset "gaugefix = :pin identifies the phase" begin
        _, pinned = Comrade.set_array(regauge(gauge_split(), :pin), arr)
        _, loose = Comrade.set_array(gauge_split(; gauge = :none), arr)
        # a wrapped angle costs two flat coordinates, so one pin takes two off
        @test Comrade.dimension(asflat(pinned)) == Comrade.dimension(asflat(loose)) - 2

        # pinning the offset by hand instead reaches the same parameterization
        _, byhand = Comrade.set_array(
            gauge_split(; offset_refant = SEFDReference(0.0)), arr
        )
        @test Comrade.dimension(asflat(pinned)) == Comrade.dimension(asflat(byhand))

        # nothing is left over, and the pinned entry no longer varies
        names = (:gpμ, :gp)
        obs = NamedTuple{names}(map(n -> getproperty(pinned, n), names))
        @test isempty(Comrade.gauge_pins(Comrade.gauge_terms(obs, names), arr))
        x, y = rand(pinned), rand(pinned)
        held = [i for i in eachindex(parent(x.gpμ)) if parent(x.gpμ)[i] == parent(y.gpμ)[i]]
        @test Comrade.sites(x.gpμ)[held] == [:AA]

        post = VLBIPosterior(skym, regauge(gauge_split(), :pin), dvis)
        @test Comrade.dimension(asflat(post)) ==
            Comrade.dimension(asflat(VLBIPosterior(skym, gauge_split(; gauge = :none), dvis))) - 2
    end

    @testset "gaugefix takes two values" begin
        m = gauge_split()
        @test_throws "must be :error or :pin" regauge(m, :ignore)
        @test regauge(m, :pin).gaugefix === :pin
        @test m.gaugefix === :error
    end

end
