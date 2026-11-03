using CairoMakie

@testset "plots" begin
    _, vis, amp, lcamp, cphase, dcoh = load_data()

    skym = SkyModel(test_model, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 64, 64))
    post = VLBIPosterior(skym, lcamp, cphase)
    θ = prior_sample(post)
    res = Comrade.residuals(post, θ)
    @test fieldplot(res[1], :uvdist, :res) isa Makie.FigureAxisPlot
    fig = Comrade.residualplot(post, θ)
    @test fig isa Makie.Figure
    @test count(c -> c isa Makie.Axis, fig.content) == 2
    @test Comrade.residualplot(post, θ; figure = (; size = (300, 400))).scene.viewport[].widths == [300, 400]

    G = Jones((x, c) -> exp(complex(x.lg, x.gp)); lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1))),
        gp = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)))
    skyp = SkyModel(test_skymodel_polarized, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 64, 64); metadata = (; lp = 0.1))
    postc = VLBIPosterior(skyp, InstrumentModel(G), dcoh)
    figc = Comrade.residualplot(postc, prior_sample(postc))
    @test count(c -> c isa Makie.Axis, figc.content) == 4
end
