using Distributions
import TransformVariables as TV
using FiniteDifferences
using Zygote
@testset "Partially fixed" begin

    d = MvLogNormal(randn(10), rand(10))
    vinds = [1, 4, 5, 6, 11, 12, 13, 14, 18, 20]
    finds = setdiff((1:20), vinds)
    dc = Comrade.PartiallyConditionedDist(d, vinds, finds, 1.0)

    @inferred logpdf(dc, rand(dc))
    @test eltype(dc) == eltype(d)
    @test length(dc) == 20

    t = asflat(dc)
    @test t isa Comrade.PartiallyFixedTransform

    @testset "PartiallyFixedTransform" begin
        tp = TV.as((gp = t,))
        x = rand(dimension(t)) .- 0.5
        p = Comrade.transform(tp, x)
        @test TV.transform(tp, TV.inverse(tp, p)).gp ≈ p.gp
        f(x) = logpdf(dc, TV.transform(tp, x).gp)
        flj(x) = (((y, lj) = TV.transform_and_logjac(tp, x));logpdf(dc, y.gp)+lj)

        fdm = central_fdm(5, 1)
        gfdf, = grad(fdm, f, x)
        gfdlj, = grad(fdm, flj, x)

        gzf,  = Zygote.gradient(f, x)
        gzflj, = Zygote.gradient(flj, x)

        @test gzf ≈ gfdf
        @test gzflj ≈ gfdlj
    end

end
