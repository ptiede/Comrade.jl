using Distributions
using Enzyme
using VLBIFiles
using VLBIImagePriors
using LinearAlgebra
using Random

function _selfcal_test_data()
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    return extract_table(uvd, Visibilities(; time_average = VLBI.GapBasedScans()))
end

function _selfcal_test_sky()
    tm(θ, meta) = θ.f1 * stretched(Gaussian(), θ.σ1, θ.σ1) + θ.f2 * shifted(Gaussian(), θ.x, 0.0)
    prior = (
        f1 = VLBIUniform(0.5, 1.5), σ1 = VLBIUniform(μas2rad(1.0), μas2rad(100.0)),
        f2 = VLBIUniform(0.1, 0.5), x = VLBIUniform(μas2rad(10.0), μas2rad(40.0)),
    )
    return SkyModel(tm, prior, imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32))
end

@instrument function amplitude_gains()
    return @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)))
        return SingleStokesGain(complex(exp(lg)))
    end
end

@instrument selfcal = true function selfcal_gains()
    return @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)))
        return SingleStokesGain(complex(exp(lg)))
    end
end

@instrument selfcal = true function named_residual()
    return @jones begin
        selfcal_residual ~ ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)))
        return SingleStokesGain(complex(exp(selfcal_residual)))
    end
end

@instrument function phase_gains()
    return @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)))
        gp ~ ArrayPrior(IIDSitePrior(IntegSeg(), VLBIUniform(-π, π)); refant = SEFDReference(0.0))
        return SingleStokesGain(exp(complex(lg, gp)))
    end
end

_withphase(θ, φ) = (sky = θ.sky, instrument = (lg = θ.instrument.lg, gp = φ))
_wrapdiff(a, b) = abs(angle(cis(a - b)))

@testset "Self-calibrated station phase" begin
    rng = Random.Xoshiro(17)
    dvis = _selfcal_test_data()
    skym = _selfcal_test_sky()
    post = VLBIPosterior(skym, selfcal_gains(), dvis)
    plain = VLBIPosterior(skym, phase_gains(), dvis)

    @testset "the likelihood is that of the plain model at φ = φ̂ + selfcal_residual" begin
        for _ in 1:3
            θ = prior_sample(rng, post)
            φ = selfcal_phase(post, θ)
            @test Comrade.loglikelihood(post, θ) ≈ Comrade.loglikelihood(plain, _withphase(θ, φ))
        end
    end

    @testset "φ̂ is the station phase of noise-free data, and zero at the fixed entries" begin
        θtrue = prior_sample(rng, plain)
        dsim = only(simulate_observation(rng, plain, θtrue; add_thermal_noise = false))
        spost = VLBIPosterior(skym, selfcal_gains(), dsim)
        θ = prior_sample(rng, spost)
        res0 = map(zero, θ.instrument.selfcal_residual)
        θ = (sky = θtrue.sky, instrument = (lg = θtrue.instrument.lg, selfcal_residual = res0))
        φ = selfcal_phase(spost, θ)
        # the true phases are referenced to a different station, so compare within each group
        d = parent(φ) .- parent(θtrue.instrument.gp)
        ref = Comrade.instrumentmodel(spost).layout.refentry
        @test maximum(_wrapdiff.(d, d[ref])) < 1.0e-8
        fixed = Comrade.fixed_indices(spost.prior.instrument.selfcal_residual)
        @test !isempty(fixed) && all(iszero, parent(φ)[fixed])
        splain = VLBIPosterior(skym, phase_gains(), dsim)
        @test Comrade.loglikelihood(spost, θ) ≈ Comrade.loglikelihood(splain, θtrue)
    end

    @testset "the gradient matches finite differences" begin
        tpost = transport_to(post, Comrade.PT.StdNormal())
        x = randn(rng, dimension(tpost)) .* 0.3
        ℓ, g = Comrade.LogDensityProblems.logdensity_and_gradient(tpost, x)
        @test isfinite(ℓ) && all(isfinite, g)
        v = normalize(randn(rng, length(x)))
        h = 1.0e-5
        fd = (Comrade.logdensityof(tpost, x .+ h .* v) - Comrade.logdensityof(tpost, x .- h .* v)) / 2h
        @test dot(g, v) ≈ fd rtol = 1.0e-4
    end

    @testset "construction is validated" begin
        @test selfcal_gains().selfcal && !amplitude_gains().selfcal
        @test_throws "the prior already has one" named_residual()
        @test_throws "needs the data" Comrade.set_array(selfcal_gains(), arrayconfig(dvis))
        dcp = extract_table(
            VLBIFiles.load(VLBIFiles.UVData, joinpath(@__DIR__, "..", "test_data.uvfits")),
            ClosurePhases()
        )
        @test_throws "exactly one table" VLBIPosterior(skym, selfcal_gains(), dcp)
    end
end
