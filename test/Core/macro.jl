using Distributions
using Random
using VLBIImagePriors

@testset "@sky macro" begin

    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 32, 32)

    @testset "parity with hand-written SkyModel" begin
        # Macro form
        @sky function macro_gauss(grid; flo = μas2rad(1.0), fhi = μas2rad(40.0))
            σ ~ VLBIUniform(flo, fhi)
            τ ~ VLBIUniform(0.35, 0.65)
            return stretched(Gaussian(), σ * τ, σ)
        end

        # Hand-written equivalent
        manual_f(θ, _, _) = stretched(Gaussian(), θ.σ * θ.τ, θ.σ)
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
        @test m_macro.f(x, m_macro.grid, m_macro.metadata) == m_manual.f(x, m_manual.grid, m_manual.metadata)
    end

    @testset "metadata flow; the grid is not metadata" begin
        @sky function with_grid(grid; ftot = 1.0)
            σ ~ VLBIUniform(μas2rad(1.0), μas2rad(40.0))
            return ftot * stretched(Gaussian(), σ, σ)
        end
        m = with_grid(g; ftot = 2.5)
        @test m.metadata.ftot == 2.5
        @test !haskey(m.metadata, :grid)
        @test m.grid === g

        x = (σ = μas2rad(10.0),)
        @test m.f(x, m.grid, m.metadata) == 2.5 * stretched(Gaussian(), x.σ, x.σ)
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
        @test m.f(x, m.grid, m.metadata) isa ContinuousImage
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
# the built-in `Jones` — it identifies param_maps by arity, never by the
# constructor's name — so a zero-arg param_map on this type is lifted and destructured too,
# and the block's priors arrive as keywords.
struct _MacroTestJones{F, P <: NamedTuple} <: Comrade.AbstractJonesMatrix
    param_map::F
    priors::P
end
_MacroTestJones(f; priors...) = _MacroTestJones(f, NamedTuple(priors))
Comrade._priors(J::_MacroTestJones) = J.priors

@testset "@instrument macro" begin
    _, dvis, _, _, _, dcoh = load_data()

    @testset "parity with hand-written InstrumentModel (StokesI)" begin
        # New syntax: one `@jones` block; the body builds and returns the Jones term, with its
        # priors as `~` lines. The natural `Ctor(entries)` form is used.
        @instrument function macro_stokesi(; refbasis = CirBasis(), gpstd = inv(π^2))
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                gp ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gpstd)); refant = SensitiveReference(0.0))
                return Jones(exp(lg + 1im * gp))
            end
        end

        G = Jones(
            (x, c) -> exp(x.lg + 1im * x.gp);
            lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
        )

        m_macro = macro_stokesi()
        m_manual = InstrumentModel(G)

        @test m_macro isa InstrumentModel
        @test keys(m_macro.prior) == keys(m_manual.prior)
        for k in keys(m_macro.prior)
            @test typeof(m_macro.prior[k]) == typeof(m_manual.prior[k])
        end

        # The synthesized param_map must be a named (serializable) function, not a closure
        @test !occursin("#", string(nameof(m_macro.jones.param_map)))

        # Same forward model on a shared parameter draw
        t = dvis
        vis = parent(Comrade.measurement(t))
        ointm, printm = bind_instrument(m_macro, t)
        ointman, _ = bind_instrument(m_manual, t)
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
                return Jones(exp(lg))
            end
        end
        @test with_kw().refbasis isa CirBasis
        @test with_kw(; refbasis = LinBasis()).refbasis isa LinBasis
        # kwarg without a default is required
        @instrument function needs_kw(; gstd)
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gstd)))
                return Jones(exp(lg))
            end
        end
        @test_throws UndefKeywordError needs_kw()
        @test needs_kw(; gstd = 0.2) isa InstrumentModel
    end

    @testset "parity with hand-written InstrumentModel (Coherencies)" begin
        @instrument function macro_pol(; refbasis = CirBasis())
            G = @jones begin
                lgR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                gpR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SensitiveReference(0.0))
                lgrat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                gprat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)), refant = SingleReference(:AA, 0.0))
                gR = exp(lgR + 1im * gpR)
                gL = gR * exp(lgrat + 1im * gprat)
                return Jones(gains(gR, gL))
            end
            D = @jones begin
                dRx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dRy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dLx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dLy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                return Jones(leakage(complex(dRx, dRy), complex(dLx, dLy)))
            end
            return G * D
        end

        m_macro = macro_pol()
        @test m_macro isa InstrumentModel
        # Priors are flat (v1) and ordered by block, then by `~` line within each block.
        @test keys(m_macro.prior) == (:lgR, :gpR, :lgrat, :gprat, :dRx, :dRy, :dLx, :dLy)

        # Both param_maps must be named (serializable)
        @test !occursin("#", string(nameof(m_macro.jones.matrices[1].param_map)))
        @test !occursin("#", string(nameof(m_macro.jones.matrices[2].param_map)))

        t = dcoh
        vis = StokesParams.(CoherencyMatrix.(collect(parent(Comrade.measurement(t))), Ref(CirBasis())))
        ointm, printm = bind_instrument(m_macro, t)
        pintm, _ = bind_instrument(InstrumentModel(), t)

        x = rand(printm)
        x.lgR .= 0; x.gpR .= 0; x.lgrat .= 0; x.gprat .= 0
        x.dRx .= 0; x.dRy .= 0; x.dLx .= 0; x.dLy .= 0
        vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
        vper = Comrade.apply_instrument(vis, pintm, (; instrument = NamedTuple()))
        @test vout ≈ vper
    end

    @testset "PerStation: priors of override-only blocks" begin
        @instrument function macro_qwp()
            G = @jones begin
                lgR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                gpR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)); refant = SingleReference(:AA, 0.0))
                lgL ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                gpL ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)); refant = SingleReference(:AA, 0.0))
                return Jones(gains(exp(complex(lgR, gpR)), exp(complex(lgL, gpL))))
            end
            Q = @jones begin
                δ ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.3)))
                return Jones(SMatrix{2, 2}(1.0 + 0im, 0, 0, cis(δ)))
            end
            return PerStation(G; AA = G * Q)
        end
        m = macro_qwp()
        @test keys(m.jones.default.priors) == (:lgR, :gpR, :lgL, :gpL)
        @test keys(m.jones.overrides.AA.matrices[2].priors) == (:δ,)
        @test keys(m.prior) == (:lgR, :gpR, :lgL, :gpL, :δ)
        _, pr = bind_instrument(m, dcoh)
        x = rand(pr)
        @test unique(sites(x.δ)) == [:AA]
        @test length(unique(sites(x.lgR))) > 1
    end

    @testset "zero-tilde parameter-free model" begin
        @instrument function paramfree(; refbasis = CirBasis())
            return @jones begin
                return Jones((x, c) -> gains(one(ComplexF64), one(ComplexF64)))
            end
        end
        m = paramfree()
        @test m isa InstrumentModel
        @test keys(m.prior) == ()
        ointm, printm = bind_instrument(m, dcoh)
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
        @test m.prior.lg === m.jones.priors.lg
        # The param_map was lifted to a named (serializable) function, not left a closure
        @test !occursin("#", string(nameof(m.jones.param_map)))
    end

    @testset "kwarg capture: param_map stays named and still captures" begin
        t = dvis
        vis = parent(Comrade.measurement(t))
        # The param_map is lifted to a *named inner function* of `_<name>_jones`, so it can
        # reference a construction-time kwarg (`off`) and remain a named, serializable function.
        @instrument function capture_kw(; refbasis = CirBasis(), off = 0.0)
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return Jones(exp(lg + off))
            end
        end
        mk = capture_kw(; off = 0.0)
        @test mk isa InstrumentModel
        @test !occursin("#", string(nameof(mk.jones.param_map)))  # named, not a gensym closure
        ointk, printk = bind_instrument(mk, t)
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
                return Jones(exp(lg) * first(map(lg -> lg, [1.0])))
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
                return Jones(round(exp(digits); digits = 2))
            end
        end
        @test kwkey() isa InstrumentModel
    end

    @testset "anonymous function in the composition keeps its scope" begin
        # A lambda in the composition that closes over a local Jones term must stay inside the
        # generated `_<name>_jones`, where that local is in scope.
        @instrument function nested_combiner(; refbasis = CirBasis())
            G = @jones begin
                lgR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return Jones(gains(exp(lgR), exp(lgR)))
            end
            D = @jones begin
                dRx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dRy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dLx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                dLy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
                return Jones(leakage(complex(dRx, dRy), complex(dLx, dLy)))
            end
            return first(map(x -> G * x, (D,)))
        end
        m = nested_combiner()
        @test m isa InstrumentModel
        @test m.jones isa Comrade.JonesProduct
        @test length(m.jones.matrices) == 2
    end

    @testset "`~` as an operator inside a body is not a misplaced prior" begin
        # `~` used as an ordinary (unary) operator in a body statement must not be misdetected as
        # a stray prior tilde — only the binary `lhs ~ rhs` form is a prior.
        @instrument function tilde_op(; refbasis = CirBasis())
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                bits = ~0x0f                       # bitwise NOT, not a prior
                scale = float(bits & 0x01)
                return Jones(exp(lg) * (1.0 + 0.0 * scale))
            end
        end
        @test tilde_op() isa InstrumentModel

        # Same for `@sky`: a body statement using `~` must be accepted.
        g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 32, 32)
        @sky function tilde_sky(grid)
            σ ~ VLBIUniform(μas2rad(1.0), μas2rad(40.0))
            bits = ~0x0f
            return (1.0 + 0.0 * (bits & 0x01)) * stretched(Gaussian(), σ, σ)
        end
        @test tilde_sky(g) isa SkyModel
    end

    @testset "standalone @jones builds a parameter-free Jones term" begin
        # A param-free `@jones` block is a valid standalone way to build a named Jones term.
        J = @jones begin
            return Jones((x, c) -> gains(one(ComplexF64), one(ComplexF64)))
        end
        @test J isa Jones
        # A standalone `@jones` with `~` priors is an error (priors have nowhere to go).
        @test_throws LoadError @eval @jones begin
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            return Jones(exp(lg))
        end
    end

    @testset "macro-expansion errors" begin
        # Positional argument on the instrument is rejected
        @test_throws LoadError @eval @instrument function haspos(array)
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return Jones(exp(lg))
            end
        end

        # `~` at the top level of the instrument body (outside a `@jones` block) is rejected
        @test_throws LoadError @eval @instrument function toplevel_tilde(; refbasis = CirBasis())
            lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
            return @jones begin
                return Jones((x, c) -> gains(one(ComplexF64), one(ComplexF64)))
            end
        end

        # `~` nested inside another expression within a `@jones` block is rejected
        @test_throws LoadError @eval @instrument function nested(; refbasis = CirBasis())
            return @jones begin
                if true
                    lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                end
                return Jones(exp(lg))
            end
        end

        # Duplicate tilde names within a block
        @test_throws LoadError @eval @instrument function dupes(; refbasis = CirBasis())
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.2)))
                return Jones(exp(lg))
            end
        end

        # Duplicate parameter name across two `@jones` blocks
        @test_throws LoadError @eval @instrument function dup_across(; refbasis = CirBasis())
            G = @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return Jones(gains(exp(lg), exp(lg)))
            end
            D = @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.2)))
                return Jones(leakage(exp(lg), exp(lg)))
            end
            return G * D
        end

        # Name clash between sampled param and kwarg
        @test_throws LoadError @eval @instrument function clash(; refbasis = CirBasis(), lg = 1.0)
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return Jones(exp(lg))
            end
        end

        # A parameterized `@jones` block whose final statement is not a single-positional-arg
        # constructor call (here a bare tuple — the `Jones` wrapper was forgotten) is rejected.
        @test_throws LoadError @eval @instrument function not_a_ctor(; refbasis = CirBasis())
            return @jones begin
                lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
                return (exp(lg), exp(lg))
            end
        end
    end
end

# --- param_map hoisting -------------------------------------------------------------------
# `@instrument` emits every param_map / combination function as a *named* function, and hoists
# it to top level when it closes over nothing. That matters for serialization: a definition
# nested inside `_<name>_jones` is a closure, so its type is gensym'd — `var"#_f_pmap_1#37"` —
# and the trailing counter is a module-global compile-time number rather than anything written
# in the source. Since the param_map is a type parameter of the Jones term, editing unrelated
# code earlier in the module shifts that counter and every previously serialized model then
# fails to load with `UndefVarError: #_f_pmap_1#37`. Hoisted, the type is the stable singleton
# `typeof(_f_pmap_1)`. A param_map that genuinely captures a construction-time keyword must
# still stay nested, or the captured value would go out of scope.
#
# These definitions are deliberately at file top level: `include` evaluates into module global
# scope, which is what allows hoisting in the first place. Defining them inside a `@testset`
# (a local scope) would leave them nested no matter what the macro does.

# Every generic function's type is named "#f", so a leading "#" distinguishes nothing; the
# closure marker is the trailing "#<counter>".
_is_hoisted(f) = isnothing(match(r"#\d+$", String(nameof(typeof(f)))))

@instrument function hoist_plain(; refbasis = CirBasis())
    return @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        gp ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        g = exp(lg + 1im * gp)
        return Jones(gains(g, g))
    end
end

# The param_map reads `gscale`, a construction-time keyword, so it cannot be hoisted.
@instrument function hoist_capture(; refbasis = CirBasis(), gscale = 2.0)
    return @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        g = exp(gscale * lg)
        return Jones(gains(g, g))
    end
end

@instrument function hoist_sandwich(; refbasis = CirBasis())
    G = @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        return Jones(gains(exp(lg), exp(lg)))
    end
    D = @jones begin
        dr ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1)))
        return Jones(leakage(complex(dr), complex(dr)))
    end
    return G * D
end

@testset "@instrument param_map hoisting" begin
    # Pin down what the marker actually distinguishes, so nothing below can pass vacuously.
    @testset "hoisted/closure probe" begin
        @test _is_hoisted(identity)
        @test !_is_hoisted(
            let y = 1
                z -> z + y
            end
        )
    end

    @testset "capture-free param_map is hoisted to top level" begin
        pm = hoist_plain().jones.param_map
        @test _is_hoisted(pm)
        @test isdefined(@__MODULE__, :_hoist_plain_pmap_1)
        @test pm === _hoist_plain_pmap_1
        @test pm((; lg = 0.3, gp = 0.2), (;)) == let g = exp(0.3 + 1im * 0.2)
            gains(g, g)
        end
    end

    @testset "capturing param_map stays nested" begin
        pm = hoist_capture().jones.param_map
        @test !_is_hoisted(pm)
        @test pm((; lg = 0.5), (;)) == let g = exp(2.0 * 0.5)
            gains(g, g)
        end
        # the capture is live: a different keyword must give a different result
        pm3 = hoist_capture(; gscale = 3.0).jones.param_map
        @test pm3((; lg = 0.5), (;)) == let g = exp(3.0 * 0.5)
            gains(g, g)
        end
    end

    @testset "param_map inside a composition is hoisted" begin
        pm = hoist_sandwich().jones.matrices[1].param_map
        @test _is_hoisted(pm)
        @test pm === _hoist_sandwich_pmap_1
    end

    @testset "hoisted type survives serialization" begin
        m = hoist_plain()
        f = tempname()
        try
            serialize(f, m.jones)
            back = deserialize(f)
            @test typeof(back) === typeof(m.jones)
            @test back.param_map === m.jones.param_map
        finally
            rm(f; force = true)
        end
    end

    # Under-collecting an assigned name is unsafe (a definition closing over it would be
    # wrongly hoisted), so every binding form must be seen — including the exception name
    # of `catch e`, which the Expr stores as a bare Symbol outside any assignment head.
    @testset "catch binding is collected as assigned" begin
        out = Symbol[]
        Comrade._collect_assigned!(
            out, :(
                try
                    f()
                catch e
                    handle(e)
                end
            )
        )
        @test :e in out
        out2 = Symbol[]
        Comrade._collect_assigned!(
            out2, :(
                try
                    f()
                catch
                end
            )
        )
        @test isempty(out2)
    end
end
