using Statistics
using Serialization
using StructArrays

@testset "PosteriorSamples" begin

    @testset "rmap over parameterless groups" begin
        # An empty group reaches `rmap` as a `Vector{@NamedTuple{}}`, not a StructArray, and
        # must pass through unchanged.
        chain = [(; sky = (; ftot = 1.0 * i, gauss = NamedTuple())) for i in 1:5]
        ps = PosteriorSamples(chain, nothing)

        m = Comrade.rmap(mean, ps)
        @test m.sky.gauss === NamedTuple()
        @test m.sky.ftot ≈ 3.0            # the real parameters still reduce

        s = Comrade.rmap(std, ps)
        @test s.sky.gauss === NamedTuple()
        @test s.sky.ftot ≈ std(1.0:5.0)

        # an empty group at the top level
        flat = [(; a = 1.0 * i, empt = NamedTuple()) for i in 1:4]
        mf = Comrade.rmap(mean, PosteriorSamples(flat, nothing))
        @test mf.empt === NamedTuple()
        @test mf.a ≈ 2.5
    end

    @testset "show with a parameterless group" begin
        # `show` builds its Mean / Std. Dev. tables through `rmap`
        chain = [(; sky = (; ftot = 1.0 * i, gauss = NamedTuple())) for i in 1:5]
        ps = PosteriorSamples(chain, nothing)
        str = sprint(show, MIME"text/plain"(), ps)
        @test occursin("PosteriorSamples", str)
        @test occursin("Mean", str)
        @test occursin("Std. Dev.", str)

        # a chain without empty groups
        plain = [(; a = 1.0 * i, b = 2.0 * i) for i in 1:5]
        @test occursin("Mean", sprint(show, MIME"text/plain"(), PosteriorSamples(plain, nothing)))
    end

    @testset "constructor call forms" begin
        chain = [(; a = 1.0 * i, b = (2.0 * i, 3.0)) for i in 1:4]
        stats = [(; lp = -1.0 * i) for i in 1:4]
        md = Dict(:sampler => :test)
        ps = PosteriorSamples(chain, stats; metadata = md)
        @test samplerstats(ps) isa StructArray
        @test samplerinfo(ps) === md
        @test PosteriorSamples(postsamples(ps), samplerstats(ps); metadata = md) == ps
        @test samplerstats(PosteriorSamples(chain, (; lp = -1.0 .* (1:4)))) == samplerstats(ps)
        @test samplerinfo(PosteriorSamples(chain, nothing)) == Dict(:sampler => :unknown)
        @test typeof(ps)(postsamples(ps), samplerstats(ps), md) == ps
        @test_deprecated PosteriorSamples(chain, stats, md)
        @test_throws "chain and stats must have the same length, got 4 and 3" PosteriorSamples(chain, stats[1:3])
    end

    @testset "resample_equal" begin
        chain = [(; a = 1.0 * i) for i in 1:4]
        ps = PosteriorSamples(chain, (; weights = [0.0, 1.0, 0.0, 0.0]))
        eq = resample_equal(ps, Int32(6))
        @test length(eq) == 6
        @test all(==(2.0), eq.a)
        @test_throws "Weights not in chain stats" resample_equal(PosteriorSamples(chain, nothing), 6)
    end

    @testset "load_samples index forms" begin
        dir = mktempdir()
        mkpath(joinpath(dir, "samples"))
        nsamples, stride = 10, 4
        chain = [(; a = 1.0 * i) for i in 1:nsamples]
        stats = [(; lp = -1.0 * i) for i in 1:nsamples]
        for (k, r) in enumerate(Iterators.partition(1:nsamples, stride))
            file = joinpath(dir, "samples", "output_scan_" * lpad(k, 5, '0') * ".jls")
            serialize(file, (samples = StructArray(chain[r]), stats = StructArray(stats[r])))
        end
        out = Comrade.DiskOutput(dir, 3, stride, nsamples)
        full = load_samples(out)
        for inds in (2:7, 3:2:9, [9, 1, 5], 4:4, isodd.(1:nsamples))
            @test postsamples(load_samples(out, inds)) == postsamples(full)[inds]
            @test load_samples(out, inds; table = :stats) == samplerstats(full)[inds]
        end
        @test_throws BoundsError load_samples(out, 0:3)
        @test_throws BoundsError load_samples(out, [11])
        @test_throws "indices must select a collection of samples" load_samples(out, 3)
    end

    @testset "chains as columns" begin
        function draws(n, μ)
            return PosteriorSamples(
                [(; sky = (; σ = μ + i, a = fill(Float64(i), 2, 2)), ρ = (1.0, 2.0)) for i in 1:n],
                [(; logp = -Float64(i)) for i in 1:n]; metadata = Dict(:sampler => :test)
            )
        end
        c1, c2 = draws(4, 0.0), draws(4, 10.0)
        ps = hcat(c1, c2)
        @test size(ps) == (4, 2)
        @test ps.sky.σ == hcat(1.0:4.0, 11.0:14.0)
        @test samplerstats(ps).logp == repeat(-(1.0:4.0), 1, 2)
        @test samplerinfo(ps) === samplerinfo(c1)
        @test ps[3, 2].sky.σ == 13.0
        @test ps[CartesianIndex(3, 2)].sky.σ == 13.0
        @test ps[7].sky.σ == 13.0
        @test ps[3, 2].sky.a == fill(3.0, 2, 2)

        second = ps[:, 2]
        @test second isa PosteriorSamples && size(second) == (4,)
        @test second.sky.σ == 11.0:14.0
        @test samplerstats(second).logp == -(1.0:4.0)
        thinned = ps[1:2:end, :]
        @test size(thinned) == (2, 2) && samplerstats(thinned).logp == repeat([-1.0, -3.0], 1, 2)

        @test reduce(hcat, [c1, c2]).sky.σ == ps.sky.σ
        @test samplerstats(reduce(hcat, [c1, c2])).logp == samplerstats(ps).logp
        joined = vcat(c1, c2)
        @test size(joined) == (8,) && samplerstats(joined).logp == [-(1.0:4.0); -(1.0:4.0)]
        @test reduce(vcat, [c1, c2]).sky.σ == joined.sky.σ

        @test Comrade.rmap(mean, ps).sky.σ ≈ 7.5
        @test vec(Comrade.rmap(x -> mean(x; dims = 1), ps).sky.σ) ≈ [2.5, 12.5]
        @test count(Returns(true), ps) == 8
        @test occursin("4 draws × 2 chains", sprint(show, MIME"text/plain"(), ps))

        @test_throws DimensionMismatch hcat(c1, draws(3, 0.0))
        @test_throws "with and without sampler statistics" hcat(c1, PosteriorSamples(postsamples(c2), nothing))
        bare = hcat(PosteriorSamples(postsamples(c1), nothing), PosteriorSamples(postsamples(c2), nothing))
        @test isnothing(samplerstats(bare)) && isnothing(samplerstats(bare[:, 1]))
    end

    @testset "load_samples reads only the files of the selected draws" begin
        dir = mktempdir()
        mkpath(joinpath(dir, "samples"))
        n, stride = 23, 10
        draws = [(; x = Float64(i), v = fill(Float64(i), 2)) for i in 1:n]
        for (f, r) in enumerate(Iterators.partition(1:n, stride))
            serialize(
                joinpath(dir, "samples", "output_scan_$(lpad(f, 5, '0')).jls"),
                (samples = Comrade.StructArray(draws[r]), stats = Comrade.StructArray([(; logp = -Float64(i)) for i in r]))
            )
        end
        serialize(joinpath(dir, "parameters.jls"), (; params = Comrade.DiskOutput(dir, 3, stride, n)))

        full = load_samples(dir)
        @test full.x == 1.0:n && samplerstats(full).logp == -(1.0:n)
        for r in (3:7, 1:4:23, 12:12, 20:-3:2)
            c = load_samples(dir, r)
            @test c.x == Float64.(r) && samplerstats(c).logp == -Float64.(r)
            @test c.v == [fill(Float64(i), 2) for i in r]
        end
        @test isnothing(samplerstats(load_samples(dir, 3:7; table = :samples)))
        @test load_samples(dir, 3:7; table = :stats).logp == -(3.0:7.0)

        @test_throws BoundsError load_samples(dir, 0:3)
        @test_throws BoundsError load_samples(dir, 20:24)
        @test_throws "no draws selected" load_samples(dir, 5:4)
        @test_throws "table must be" load_samples(dir; table = :chain)

        write(joinpath(dir, "samples", "output_scan_00002.jls"), "not a sample file")
        @test load_samples(dir, 1:20:21).x == [1.0, 21.0]
    end

end
