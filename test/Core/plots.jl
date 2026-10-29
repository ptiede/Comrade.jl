using CairoMakie

@testset "plots" begin
    _, vis, amp, lcamp, cphase, dcoh = load_data()
    rows = datarows(vis)
    @test length(rows) == count(>(0), parent(weight(vis)))
    @test rows.uvdist ≈ hypot.(rows.U, rows.V)
    crow = datarows(cphase)
    @test length(crow) == length(cphase)
    @test all(>(0), crow.uvdist)
    @test datarows(dcoh).measurement isa AbstractVector{<:AbstractMatrix}

    for t in (vis, amp, lcamp, cphase)
        @test plotfields(t, :uvdist, :measurement) isa Makie.FigureAxisPlot
    end
    fig = plotfields(dcoh, :uvdist, :measurement)
    @test fig isa Makie.Figure
    @test count(c -> c isa Makie.Axis, fig.content) == 4
    @test plotfields(dcoh, (:AA, :LM), :Ti, :amp) isa Makie.Figure
    @test plotfields(vis, :U, :V) isa Makie.FigureAxisPlot
    @test plotfields(vis, (:AA, :LM), :Ti, :amp) isa Makie.FigureAxisPlot
    @test plotfields(dcoh, :uvdist, x -> abs(x.measurement[1, 1])) isa Makie.FigureAxisPlot
    @test_throws "closures are not on one baseline" plotfields(cphase, (:AA, :LM), :Ti, :measurement)

    fig = Figure()
    ax, = baselineplot(fig[1, 1], cphase, :uvdist, :measwnoise)
    @test baselineplot!(ax, vis, (:AA, :LM), :uvdist, :amp) isa Makie.AbstractPlot

    skym = SkyModel(test_model, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 64, 64))
    post = VLBIPosterior(skym, lcamp, cphase)
    θ = prior_sample(post)
    res = Comrade.residuals(post, θ)
    @test plotfields(res[1], :uvdist, :res) isa Makie.FigureAxisPlot
    fig = Comrade.residual(post, θ)
    @test fig isa Makie.Figure
    @test count(c -> c isa Makie.Axis, fig.content) == 2

    G = Jones((x, c) -> exp(complex(x.lg, x.gp)); lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1))),
        gp = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)))
    skyp = SkyModel(test_skymodel_polarized, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 64, 64); metadata = (; lp = 0.1))
    postc = VLBIPosterior(skyp, InstrumentModel(G), dcoh)
    figc = Comrade.residual(postc, prior_sample(postc))
    @test count(c -> c isa Makie.Axis, figc.content) == 4
end
