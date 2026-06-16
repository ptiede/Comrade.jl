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
        @test m.prior.ρs isa NTuple{3, <:VLBIImagePriors.PushforwardDistribution}
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

# A user-defined param-bearing Jones type. The `@instrument` macro must treat it exactly like
# the built-ins (`JonesG`, `JonesD`, ...) — it identifies param_maps by arity, never by the
# constructor's name — so a zero-arg param_map on this type is lifted and destructured too.
struct _MacroTestJones{F} <: Comrade.AbstractJonesMatrix
    param_map::F
end

@testset "@instrument macro" begin
    _, dvis, _, _, _, dcoh = load_data()

    @testset "parity with hand-written InstrumentModel (StokesI)" begin
        # New syntax: one `@jones` block; the body builds and returns the Jones term, with its
        # priors as `~` lines. The natural `Ctor(entries)` form is used.
        @instrument function macro_stokesi(; refbasis = CirBasis(), gpstd = inv(π^2))
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                gp ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gpstd)); refant = SEFDReference(0.0))
                return SingleStokesGain(exp(lg + 1im * gp))
            end
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
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gstd)))
                return SingleStokesGain(exp(lg))
            end
        end
        @test with_kw().refbasis isa CirBasis
        @test with_kw(; refbasis = LinBasis()).refbasis isa LinBasis
        # kwarg without a default is required
        @instrument function needs_kw(; gstd)
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gstd)))
                return SingleStokesGain(exp(lg))
            end
        end
        @test_throws UndefKeywordError needs_kw()
        @test needs_kw(; gstd = 0.2) isa InstrumentModel
    end

    @testset "parity with hand-written InstrumentModel (Coherencies)" begin
        @instrument function macro_pol(; refbasis = CirBasis())
            G = @jones begin
                lgR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                gpR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); phase = true, refant = SEFDReference(0.0))
                lgrat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)), phase = false)
                gprat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)), refant = SingleReference(:AA, 0.0))
                gR = exp(lgR + 1im * gpR)
                gL = gR * exp(lgrat + 1im * gprat)
                return JonesG((gR, gL))
            end
            D = @jones begin
                dRx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dRy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dLx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dLy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                return JonesD((complex(dRx, dRy), complex(dLx, dLy)))
            end
            R = @jones begin
                return JonesR(; add_fr = true)
            end
            return JonesSandwich(G, D, R) do g, d, r
                return g * d * r
            end
        end

        m_macro = macro_pol()
        @test m_macro isa InstrumentModel
        # Priors are flat (v1) and ordered by block, then by `~` line within each block.
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
            return @jones begin
                return JonesR(; add_fr = true)
            end
        end
        m = responseonly()
        @test m isa InstrumentModel
        @test keys(m.prior) == ()
        ointm, printm = Comrade.set_array(m, arrayconfig(dcoh))
        @test ointm isa Comrade.ObservedInstrumentModel
    end

    @testset "name-agnostic: user-defined param-bearing Jones type" begin
        # No constructor-name list: the macro identifies the param_map structurally (the single
        # positional argument of the returned constructor), so a user-defined `AbstractJonesMatrix`
        # is handled exactly like the built-ins.
        @instrument function customjones(; refbasis = CirBasis())
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return _MacroTestJones(exp(lg))
            end
        end
        m = customjones()
        @test m isa InstrumentModel
        @test m.jones isa _MacroTestJones
        @test keys(m.prior) == (:lg,)
        # The param_map was lifted to a named (serializable) function, not left a closure
        @test !occursin("#", string(nameof(m.jones.param_map)))
    end

    @testset "kwarg capture: param_map stays named and still captures" begin
        vis = Comrade.measurement(dvis)
        # The param_map is lifted to a *named inner function* of `_<name>_jones`, so it can
        # reference a construction-time kwarg (`off`) and remain a named, serializable function.
        @instrument function capture_kw(; refbasis = CirBasis(), off = 0.0)
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return SingleStokesGain(exp(lg + off))
            end
        end
        mk = capture_kw(; off = 0.0)
        @test mk isa InstrumentModel
        @test !occursin("#", string(nameof(mk.jones.param_map)))  # named, not a gensym closure
        ointk, printk = Comrade.set_array(mk, arrayconfig(dvis))
        xk = rand(printk)
        xk.lg .= 0
        @test Comrade.apply_instrument(vis, ointk, (; instrument = xk)) ≈ vis
    end

    @testset "in-param_map shadowing is Julia's job" begin
        # Inside a `@jones` block ordinary Julia scoping applies, so a nested closure may freely
        # reuse a sampled-parameter name: the outer `lg` is the destructured param and the
        # inner `lg -> lg` is the nested closure's argument.
        @instrument function nested_closure(; refbasis = CirBasis())
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return SingleStokesGain(exp(lg) * first(map(lg -> lg, [1.0])))
            end
        end
        mn = nested_closure()
        @test mn isa InstrumentModel
        @test !occursin("#", string(nameof(mn.jones.param_map)))  # lifted to a named function

        # A keyword-argument KEY inside a param_map that matches a parameter name is fine: the
        # param `digits` is destructured and used in `exp(digits)`, while the `digits = 2`
        # keyword key is plain Julia syntax, untouched.
        @instrument function kwkey(; refbasis = CirBasis())
            return @jones begin
                digits ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return SingleStokesGain(round(exp(digits); digits = 2))
            end
        end
        @test kwkey() isa InstrumentModel
    end

    @testset "nested anonymous function in the combination function keeps its scope" begin
        # A lambda nested inside the lifted `JonesSandwich` combination function may close over
        # the combinator's own arguments. The lifted (named) outer function must keep that nested
        # lambda inline; hoisting it to a sibling named function would drop `g`/`d` from scope and
        # crash at evaluation time.
        @instrument function nested_combiner(; refbasis = CirBasis())
            G = @jones begin
                lgR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return JonesG((exp(lgR), exp(lgR)))
            end
            D = @jones begin
                dRx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dRy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dLx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dLy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                return JonesD((complex(dRx, dRy), complex(dLx, dLy)))
            end
            R = @jones begin
                return JonesR(; add_fr = true)
            end
            # The `map(x -> g * x, ...)` lambda closes over `g`, an argument of the do-block.
            return JonesSandwich(G, D, R) do g, d, r
                return first(map(x -> g * x, (d,))) * r
            end
        end
        m = nested_combiner()
        @test m isa InstrumentModel
        # Combination function is still named (serializable)
        @test !occursin("#", string(nameof(m.jones.jones_map)))
        # Directly exercise the lifted combination function: the nested `x -> g * x` lambda
        # captures `g`, the do-block's argument. Before the fix this lambda was hoisted to a
        # sibling named function and threw UndefVarError(:g) on the first call. `JonesSandwich`
        # wraps the combinator in `Base.Splat`, so it is called with a single tuple of matrices.
        jm = m.jones.jones_map
        @test jm((2.0, 3.0, 5.0)) == (2.0 * 3.0) * 5.0
    end

    @testset "`~` as an operator inside a body is not a misplaced prior" begin
        # `~` used as an ordinary (unary) operator in a body statement must not be misdetected as
        # a stray prior tilde — only the binary `lhs ~ rhs` form is a prior.
        @instrument function tilde_op(; refbasis = CirBasis())
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                bits = ~0x0f                       # bitwise NOT, not a prior
                scale = float(bits & 0x01)
                return SingleStokesGain(exp(lg) * (1.0 + 0.0 * scale))
            end
        end
        @test tilde_op() isa InstrumentModel

        # Same for `@sky`: a body statement using `~` must be accepted.
        g = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)
        @sky function tilde_sky(grid)
            σ ~ VLBIUniform(μas2rad(1.0), μas2rad(40.0))
            bits = ~0x0f
            return (1.0 + 0.0 * (bits & 0x01)) * stretched(Gaussian(), σ, σ)
        end
        @test tilde_sky(g) isa SkyModel
    end

    @testset "standalone @jones builds a parameter-free Jones term" begin
        # A param-free `@jones` block is a valid standalone way to build a named Jones term.
        R = @jones begin
            return JonesR(; add_fr = true)
        end
        @test R isa JonesR
        # A standalone `@jones` with `~` priors is an error (priors have nowhere to go).
        @test_throws LoadError @eval @jones begin
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            return SingleStokesGain(exp(lg))
        end
    end

    @testset "macro-expansion errors" begin
        # Positional argument on the instrument is rejected
        @test_throws LoadError @eval @instrument function haspos(array)
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return SingleStokesGain(exp(lg))
            end
        end

        # `~` at the top level of the instrument body (outside a `@jones` block) is rejected
        @test_throws LoadError @eval @instrument function toplevel_tilde(; refbasis = CirBasis())
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            return @jones begin
                return JonesR(; add_fr = true)
            end
        end

        # `~` nested inside another expression within a `@jones` block is rejected
        @test_throws LoadError @eval @instrument function nested(; refbasis = CirBasis())
            return @jones begin
                if true
                    lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                end
                return SingleStokesGain(exp(lg))
            end
        end

        # Duplicate tilde names within a block
        @test_throws LoadError @eval @instrument function dupes(; refbasis = CirBasis())
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.2)))
                return SingleStokesGain(exp(lg))
            end
        end

        # Duplicate parameter name across two `@jones` blocks
        @test_throws LoadError @eval @instrument function dup_across(; refbasis = CirBasis())
            G = @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return JonesG((exp(lg), exp(lg)))
            end
            D = @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.2)))
                return JonesD((exp(lg), exp(lg)))
            end
            return JonesSandwich(G, D) do g, d
                return g * d
            end
        end

        # Name clash between sampled param and kwarg
        @test_throws LoadError @eval @instrument function clash(; refbasis = CirBasis(), lg = 1.0)
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return SingleStokesGain(exp(lg))
            end
        end

        # A parameterized `@jones` block whose final statement is not a single-positional-arg
        # constructor call (here a bare tuple — the `JonesG` wrapper was forgotten) is rejected.
        @test_throws LoadError @eval @instrument function not_a_ctor(; refbasis = CirBasis())
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return (exp(lg), exp(lg))
            end
        end
    end
end
