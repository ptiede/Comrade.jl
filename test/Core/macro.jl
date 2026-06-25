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
        x = rand(rng, Comrade.HypercubeTransform.NamedDist(m_macro.prior))
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
        x = rand(Comrade.HypercubeTransform.NamedDist(m.prior))
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
        x = rand(Comrade.HypercubeTransform.NamedDist(m.prior))
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
      
@testset "@instrument macro" begin
    _, dvis, _, _, _, dcoh = load_data()

    @testset "parity with hand-written InstrumentModel (StokesI)" begin
        # New syntax: reference the sampled-parameter names directly (no dummy `x`).
        @instrument function macro_stokesi(; refbasis = CirBasis(), gpstd = inv(π^2))
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            gp ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gpstd)); refant = SEFDReference(0.0))
            return SingleStokesGain(exp(lg + 1im * gp))
        end

        G = SingleStokesGain(x -> exp(x.lg + 1im * x.gp))
        manual_prior = (
            lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SEFDReference(0.0)),
        )

        m_macro = macro_stokesi()
        m_manual = InstrumentModel(G, manual_prior)

        @test m_macro isa InstrumentModel
        @test keys(m_macro.prior) == keys(m_manual.prior)
        for k in keys(m_macro.prior)
            @test typeof(m_macro.prior[k]) == typeof(m_manual.prior[k])
        end

        # The synthesized param_map must be a named (serializable) function, not a closure
        @test !occursin("#", string(nameof(m_macro.jones.param_map)))

        # Same forward model on a shared parameter draw
        vis = Comrade.measurement(dvis)
        ointm, printm = Comrade.set_array(m_macro, arrayconfig(dvis))
        ointman, _ = Comrade.set_array(m_manual, arrayconfig(dvis))
        x = rand(printm)
        vout_macro = Comrade.apply_instrument(vis, ointm, (; instrument = x))
        vout_manual = Comrade.apply_instrument(vis, ointman, (; instrument = x))
        @test vout_macro ≈ vout_manual

        # Identity params round-trip back to the input visibilities
        x.lg .= 0
        x.gp .= 0
        @test Comrade.apply_instrument(vis, ointm, (; instrument = x)) ≈ vis
    end

    @testset "refbasis and kwarg flow" begin
        @instrument function with_kw(; refbasis = CirBasis(), gstd = 0.1)
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gstd)))
            return SingleStokesGain(exp(lg))
        end
        @test with_kw().refbasis isa CirBasis
        @test with_kw(; refbasis = LinBasis()).refbasis isa LinBasis
        # kwarg without a default is required
        @instrument function needs_kw(; gstd)
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gstd)))
            return SingleStokesGain(exp(lg))
        end
        @test_throws UndefKeywordError needs_kw()
        @test needs_kw(; gstd = 0.2) isa InstrumentModel
    end

    @testset "parity with hand-written InstrumentModel (Coherencies)" begin
        @instrument function macro_pol(; refbasis = CirBasis())
            lgR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            gpR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); phase = true, refant = SEFDReference(0.0))
            lgrat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)), phase = false)
            gprat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)), refant = SingleReference(:AA, 0.0))
            dRx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
            dRy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
            dLx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
            dLy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))

            # zero-arg do-block (multi-line) and bare expression (one-liner), no dummy `x`
            G = JonesG() do
                gR = exp(lgR + 1im * gpR)
                gL = gR * exp(lgrat + 1im * gprat)
                return gR, gL
            end
            D = JonesD((complex(dRx, dRy), complex(dLx, dLy)))
            R = JonesR(; add_fr = true)
            return JonesSandwich(G, D, R) do g, d, r
                return g * d * r
            end
        end

        m_macro = macro_pol()
        @test m_macro isa InstrumentModel
        @test keys(m_macro.prior) == (:lgR, :gpR, :lgrat, :gprat, :dRx, :dRy, :dLx, :dLy)

        # Both param_maps and the combination function must be named (serializable)
        @test !occursin("#", string(nameof(m_macro.jones.matrices[1].param_map)))
        @test !occursin("#", string(nameof(m_macro.jones.matrices[2].param_map)))
        @test !occursin("#", string(nameof(m_macro.jones.jones_map)))

        vis = CoherencyMatrix.(Comrade.measurement(dcoh), Ref(CirBasis()))
        ointm, printm = Comrade.set_array(m_macro, arrayconfig(dcoh))
        pintm, _ = Comrade.set_array(InstrumentModel(JonesR(; add_fr = true)), arrayconfig(dcoh))

        x = rand(printm)
        x.lgR .= 0; x.gpR .= 0; x.lgrat .= 0; x.gprat .= 0
        x.dRx .= 0; x.dRy .= 0; x.dLx .= 0; x.dLy .= 0
        vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
        vper = Comrade.apply_instrument(vis, pintm, (; instrument = NamedTuple()))
        @test vout ≈ vper
    end

    @testset "zero-tilde response-only model" begin
        @instrument function responseonly(; refbasis = CirBasis())
            return JonesR(; add_fr = true)
        end
        m = responseonly()
        @test m isa InstrumentModel
        @test keys(m.prior) == ()
        ointm, printm = Comrade.set_array(m, arrayconfig(dcoh))
        @test ointm isa Comrade.ObservedInstrumentModel
    end

    @testset "explicit closure escape hatch and kwarg capture" begin
        # Explicit `x -> ...` is left untouched (required for custom Jones types)
        @instrument function explicit_x(; refbasis = CirBasis())
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            return SingleStokesGain(exp(lg))
        end
        mx = explicit_x()
        @test mx isa InstrumentModel
        vis = Comrade.measurement(dvis)
        ointx, printx = Comrade.set_array(mx, arrayconfig(dvis))
        xx = rand(printx)
        xx.lg .= 0
        @test Comrade.apply_instrument(vis, ointx, (; instrument = xx)) ≈ vis

        # A param_map that captures a construction-time kwarg stays a closure but still works
        @instrument function capture_kw(; refbasis = CirBasis(), off = 0.0)
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            return SingleStokesGain(exp(lg + off))
        end
        mk = capture_kw(; off = 0.0)
        @test mk isa InstrumentModel
        ointk, printk = Comrade.set_array(mk, arrayconfig(dvis))
        xk = rand(printk)
        xk.lg .= 0
        @test Comrade.apply_instrument(vis, ointk, (; instrument = xk)) ≈ vis
    end

    @testset "macro-expansion errors" begin
        # Positional argument is rejected
        @test_throws LoadError @eval @instrument function haspos(array)
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            return SingleStokesGain(exp(lg))
        end

        # Tilde nested inside another expression is rejected
        @test_throws LoadError @eval @instrument function nested(; refbasis = CirBasis())
            if true
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            end
            return SingleStokesGain(exp(lg))
        end

        # Duplicate tilde names
        @test_throws LoadError @eval @instrument function dupes(; refbasis = CirBasis())
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.2)))
            return SingleStokesGain(exp(lg))
        end

        # Name clash between sampled param and kwarg
        @test_throws LoadError @eval @instrument function clash(; refbasis = CirBasis(), lg = 1.0)
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            return SingleStokesGain(exp(lg))
        end
    end
end
