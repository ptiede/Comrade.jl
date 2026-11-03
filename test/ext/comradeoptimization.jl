using Comrade, Optimization
using OptimizationLBFGSB, Distributions, VLBIImagePriors
using OptimizationBBO, OptimizationGCMAES, OptimizationOptimisers
using Enzyme
using Test

const COExt = Base.get_extension(Comrade, :ComradeOptimizationExt)
# The cube (`ascube`) transform is a `TransportedDistribution` over a `StdUniform` latent.
const Ahyper = Comrade.CubeTransport


@testset "ComradeOptimizationExt" begin
    _, _, _, lcamp, cphase = load_data()
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(test_model, test_prior(), g)
    post = VLBIPosterior(skym, lcamp, cphase; admode = set_runtime_activity(Enzyme.Reverse))
    tpost = asflat(post)
    x0 = transform(
        tpost, [
            0.0,
            -0.4,
            0.0,
            2.0,
            0.0,
            -1.0,
            0.0,
            0.5,
            2.0,
            2.0,
        ]
    )

    # This is a dummy test for the heuristic to ensure that for large problems we switch to asflat
    skymbig = SkyModel(test_model, merge(test_prior(), (; foo = MvNormal(ones(100)))), g)
    postbig = VLBIPosterior(skymbig, lcamp, cphase; admode = set_runtime_activity(Enzyme.Reverse))

    @testset "Heuristic Test" begin
        @test !isa(COExt._transport_heuristic(post, LBFGSB(), nothing, nothing).transform, Ahyper)
        @test isa(COExt._transport_heuristic(post, LBFGSB(), nothing, StdUniform()).transform, Ahyper)
        @test !isa(COExt._transport_heuristic(post, LBFGSB(), nothing, TVFlat()).transform, Ahyper)
        @test isa(COExt._transport_heuristic(post, BBO_adaptive_de_rand_1_bin_radiuslimited(), nothing, nothing).transform, Ahyper)
        @test isa(COExt._transport_heuristic(post, GCMAESOpt(), nothing, nothing).transform, Ahyper)
        @test_throws "You are using the" COExt._transport_heuristic(postbig, GCMAESOpt(), nothing, nothing)
        # We don't actually use the bounds when constructing so just test that the condition in the heuristic is correct
        @test !isa(COExt._transport_heuristic(postbig, GCMAESOpt(), prior_sample(postbig), nothing).transform, Ahyper)
        @test !isa(COExt._transport_heuristic(postbig, OptimizationOptimisers.Adam(), prior_sample(postbig), nothing).transform, Ahyper)
    end

    xopt2, sol = comrade_opt(post, LBFGSB(); initial_params = x0, maxiters = 10_000)

    @testset "transport_method" begin
        xc, _ = comrade_opt(post, LBFGSB(); initial_params = x0, transport_method = StdUniform(), maxiters = 10)
        @test keys(xc) == keys(xopt2)
        xd, _ = @test_deprecated comrade_opt(post, LBFGSB(); initial_params = x0, transform = ascube, maxiters = 10)
        @test keys(xd) == keys(xopt2)
        @test_throws "pass only `transport_method`" comrade_opt(post, LBFGSB(); transform = ascube, transport_method = StdUniform())
    end

    xopt = xopt2.sky
    @test isapprox(xopt.f1 / xopt.f2, 2.0, atol = 1.0e-3)
    @test isapprox(xopt.σ1 * 2 * sqrt(2 * log(2)), μas2rad(40.0), rtol = 1.0e-3)
    @test isapprox(xopt.σ1 * xopt.τ1 * 2 * sqrt(2 * log(2)), μas2rad(20.0), rtol = 1.0e-3)
    @test isapprox(xopt.ξ1, π / 3, atol = 1.0e-3)
    @test isapprox(xopt.σ2 * 2 * sqrt(2 * log(2)), μas2rad(20.0), atol = 1.0e-3)
    @test isapprox(xopt.σ2 * xopt.τ2 * 2 * sqrt(2 * log(2)), μas2rad(10.0), rtol = 1.0e-3)
    @test isapprox(xopt.ξ2, π / 6, atol = 1.0e-3)
    @test isapprox(xopt.x, μas2rad(30.0), rtol = 1.0e-3)
    @test isapprox(xopt.y, μas2rad(30.0), rtol = 1.0e-3)
    @test sum(chi2(post, xopt2)) < 0.1

end
