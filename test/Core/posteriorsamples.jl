using Statistics
using Serialization
using StructArrays

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
        for inds in (2:7, 3:2:9, [9, 1, 5], 4:4, Int[], isodd.(1:nsamples))
            @test postsamples(load_samples(out, inds)) == postsamples(full)[inds]
            @test load_samples(out, inds; table = :stats) == samplerstats(full)[inds]
        end
        @test_throws BoundsError load_samples(out, 0:3)
        @test_throws BoundsError load_samples(out, [11])
        @test_throws "indices must select a collection of samples" load_samples(out, 3)
    end

end
