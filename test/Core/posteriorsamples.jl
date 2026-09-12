using Statistics

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

end
