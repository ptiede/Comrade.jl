using Distributions
using Random
using VLBIImagePriors

@testset "@sky macro" begin

    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)

    @testset "parity with hand-written SkyModel" begin
        # Macro form
        @sky function macro_gauss(grid; flo = μas2rad(1.0), fhi = μas2rad(40.0))
            σ ~ VLBIUniform(flo, fhi)
            τ ~ VLBIUniform(0.35, 0.65)
            return stretched(Gaussian(), σ * τ, σ)
        end

        # Hand-written equivalent
        manual_f(θ, meta) = stretched(Gaussian(), θ.σ * θ.τ, θ.σ)
        manual_prior = (
            σ = VLBIUniform(μas2rad(1.0), μas2rad(40.0)),
            τ = VLBIUniform(0.35, 0.65),
        )

        m_macro = macro_gauss(g)
        m_manual = SkyModel(manual_f, manual_prior, g)

        @test m_macro isa SkyModel
        @test m_manual isa SkyModel

        # Same prior shape
        @test keys(m_macro.prior) == keys(m_manual.prior)
        for k in keys(m_macro.prior)
            @test typeof(m_macro.prior[k]) == typeof(m_manual.prior[k])
        end

        # Same model output on a fixed θ
        rng = Random.Xoshiro(42)
        x = rand(rng, Comrade.NamedDist(m_macro.prior))
        @test m_macro.f(x, m_macro.metadata) == m_manual.f(x, m_manual.metadata)
    end

    @testset "metadata flow and grid in metadata" begin
        @sky function with_grid(grid; ftot = 1.0)
            σ ~ VLBIUniform(μas2rad(1.0), μas2rad(40.0))
            return ftot * stretched(Gaussian(), σ, σ)
        end
        m = with_grid(g; ftot = 2.5)
        @test m.metadata.ftot == 2.5
        @test m.metadata.grid === g
        @test m.grid === g

        x = (σ = μas2rad(10.0),)
        @test m.f(x, m.metadata) == 2.5 * stretched(Gaussian(), x.σ, x.σ)
    end

    @testset "kwarg without default is required" begin
        @sky function needs_kw(grid; ftot)
            σ ~ VLBIUniform(μas2rad(1.0), μas2rad(40.0))
            return ftot * stretched(Gaussian(), σ, σ)
        end
        @test_throws UndefKeywordError needs_kw(g)
        m = needs_kw(g; ftot = 1.5)
        @test m.metadata.ftot == 1.5
    end

    @testset "tuple-of-IID priors via ntuple" begin
        @sky function ntuple_prior(grid; n = 3, lo = 0.01, hi = 10.0)
            ρs ~ ntuple(Returns(VLBIUniform(lo, hi)), n)
            return stretched(Gaussian(), ρs[1], ρs[2])
        end
        m = ntuple_prior(g)
        @test m.prior.ρs isa NTuple{3, <:VLBIImagePriors.AffineDistribution}
        x = rand(Comrade.NamedDist(m.prior))
        @test x.ρs isa NTuple{3, Float64}
    end

    @testset "hierarchical prior on tilde RHS" begin
        cprior = corr_image_prior(g, μas2rad(20.0))
        @test cprior isa VLBIImagePriors.HierarchicalPrior

        @sky function hier_imager(grid; cprior, σscale = 0.1)
            c ~ cprior
            σimg ~ VLBIExponential(σscale)
            rast = to_simplex(CenteredLR(), σimg .* c.params)
            return ContinuousImage(rast, grid, BSplinePulse{3}())
        end
        m = hier_imager(g; cprior = cprior)
        @test m.prior.c === cprior
        x = rand(Comrade.NamedDist(m.prior))
        @test haskey(x, :c) && haskey(x, :σimg)
        @test m.f(x, m.metadata) isa ContinuousImage
    end

    @testset "VLBIPosterior compatibility" begin
        @sky function gauss2(grid)
            σ ~ VLBIUniform(μas2rad(1.0), μas2rad(40.0))
            τ ~ VLBIUniform(0.35, 0.65)
            f1 ~ VLBIUniform(0.8, 1.2)
            return f1 * stretched(Gaussian(), σ * τ, σ)
        end
        m = gauss2(g)
        # Smoke-test through the full posterior pipeline
        _, vis, _, _, _ = load_data()
        post = VLBIPosterior(m, vis)
        x = prior_sample(post)
        @test isfinite(logdensityof(post, x))
        tpostf = asflat(post)
        tpostc = ascube(post)
        @test isfinite(logdensityof(tpostf, prior_sample(tpostf)))
        @test isfinite(logdensityof(tpostc, prior_sample(tpostc)))
    end

    @testset "macro-expansion errors" begin
        # Tilde nested inside another expression is rejected
        @test_throws LoadError @eval @sky function nested_tilde(grid)
            if true
                σ ~ VLBIUniform(0.0, 1.0)
            end
            return Gaussian()
        end

        # Name clash between sampled param and kwarg
        @test_throws LoadError @eval @sky function clash(grid; σ = 1.0)
            σ ~ VLBIUniform(0.0, 1.0)
            return stretched(Gaussian(), σ, σ)
        end

        # Name clash with positional grid name
        @test_throws LoadError @eval @sky function gridclash(grid)
            grid ~ VLBIUniform(0.0, 1.0)
            return Gaussian()
        end

        # No tildes at all
        @test_throws LoadError @eval @sky function notildes(grid)
            return Gaussian()
        end

        # Multiple positional args rejected
        @test_throws LoadError @eval @sky function twopos(grid, extra)
            σ ~ VLBIUniform(0.0, 1.0)
            return Gaussian()
        end

        # Duplicate tilde names
        @test_throws LoadError @eval @sky function dupes(grid)
            σ ~ VLBIUniform(0.0, 1.0)
            σ ~ VLBIUniform(0.0, 2.0)
            return Gaussian()
        end
    end
end
