using Statistics, Serialization

@testset "PosteriorSamples" begin

    @testset "rmap over parameterless groups" begin
        # A model component that is switched off contributes a group with no parameters, so
        # every draw carries an empty NamedTuple (e.g. a sky model's `gauss` when no extra
        # Gaussian is fitted). StructArrays cannot split that into component arrays, so it
        # reaches `rmap` as a plain `Vector{@NamedTuple{}}` rather than a `StructArray`, and
        # the generic `rmap(f, x) = f(x)` fallback would call the reduction on it —
        # `mean` then fails in `sum(...)/n` with `no method matching /(::@NamedTuple{}, ::Int64)`.
        # There is nothing to average, so the group must pass through untouched.
        chain = [(; sky = (; ftot = 1.0 * i, gauss = NamedTuple())) for i in 1:5]
        ps = PosteriorSamples(chain, nothing)

        m = Comrade.rmap(mean, ps)
        @test m.sky.gauss === NamedTuple()
        @test m.sky.ftot ≈ 3.0            # the real parameters still reduce

        s = Comrade.rmap(std, ps)
        @test s.sky.gauss === NamedTuple()
        @test s.sky.ftot ≈ std(1.0:5.0)

        # an empty group at the top level, not just nested
        flat = [(; a = 1.0 * i, empt = NamedTuple()) for i in 1:4]
        mf = Comrade.rmap(mean, PosteriorSamples(flat, nothing))
        @test mf.empt === NamedTuple()
        @test mf.a ≈ 2.5
    end

    @testset "show with a parameterless group" begin
        # `show` builds its Mean / Std. Dev. tables through `rmap`, so an empty group used to
        # take the whole display down.
        chain = [(; sky = (; ftot = 1.0 * i, gauss = NamedTuple())) for i in 1:5]
        ps = PosteriorSamples(chain, nothing)
        str = sprint(show, MIME"text/plain"(), ps)
        @test occursin("PosteriorSamples", str)
        @test occursin("Mean", str)
        @test occursin("Std. Dev.", str)

        # and it is unaffected for a chain with no empty groups
        plain = [(; a = 1.0 * i, b = 2.0 * i) for i in 1:5]
        @test occursin("Mean", sprint(show, MIME"text/plain"(), PosteriorSamples(plain, nothing)))
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
