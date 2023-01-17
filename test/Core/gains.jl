using Distributions
using Plots
using Tables

@testset "gains" begin
    _,vis, amp, lcamp, cphase = load_data()

    st = scantable(vis)
    stamp = scantable(amp)
    stcp = scantable(cphase)
    stlca = scantable(lcamp)

    gcache = GainCache(st)


    θ = (f1 = 1.0,
         σ1 = 1.2e-10,
         τ1 = 0.25,
         ξ1 = -0.42,
         f2 = 0.5,
         σ2 = 8.0e-12,
         τ2 = 0.4,
         ξ2 = 0.0,
         x = 4.0e-11,
         y = 1.5e-10
        )

    m = test_model(θ)

    tel = stations(vis)
    gamp_prior = NamedTuple{Tuple(tel)}(ntuple(_->LogNormal(0.0, 0.1), length(tel)))
    gph_prior = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))

    gamp = CalPrior(gamp_prior, st)
    gpha = CalPrior(gph_prior, st)

    ga = fill(1.0, size(rand(gamp))...)
    gm = JonesModel(jonesStokes(ga, gcache), m)

    c1 = caltable(gm)
    c2 = caltable(gcache, ga)
    c3 = caltable(amp, ga)

    @testset "caltable test" begin
        @test Tables.istable(typeof(c1))
        @test Tables.rowaccess(typeof(c1))
        @test Tables.rows(c1) === c1
        @test Tables.columnaccess(c1)
        clmns = Tables.columns(c1)
        @test clmns[1] == Comrade.scantimes(c1)
        @test Bool(prod(skipmissing(Tables.matrix(clmns)[:,begin+1:end]) .== skipmissing(Comrade.gmat(c1))))
        @test c1.time == Comrade.scantimes(c1)
        @test c1.time == Tables.getcolumn(c1, 1)
        @test c1.AA == Tables.getcolumn(c1, :AA)
        @test c1.AA == Tables.getcolumn(c1, 2)
        @test Tables.columnnames(c1) == [:time, sort(stations(amp))...]

        c1row = first(c1)
        @test eltype(c1) == typeof(c1row)
        @test c1row.time == c1.time[1]
        @test c1row.AA == c1.AA[1]
        @test Tables.getcolumn(c1row, :AA) == c1.AA[1]
        @test Tables.getcolumn(c1row, :time) == c1.time[1]
        @test Tables.getcolumn(c1row, 2) == c1.AA[1]
        @test Tables.getcolumn(c1row, 1) == c1.time[1]
        @test propertynames(c1) == propertynames(c1row) == [:time, sort(stations(amp))...]


    end


    @test prod(skipmissing(Comrade.gmat(c3) .== Comrade.gmat(c1)))
    @test prod(skipmissing(Comrade.gmat(c3) .== Comrade.gmat(c2)))
    @test prod(skipmissing(Comrade.gmat(c2) .== Comrade.gmat(c1)))

    plot(c3, layout=(3,3), size=(600,500))

    ac = arrayconfig(vis)
    cpac = arrayconfig(cphase)
    lcac = arrayconfig(lcamp)
    @test visibilities(gm, ac) ≈ visibilities(m, ac)
    @test amplitudes(gm, ac) ≈ amplitudes(m, ac)
    @test closure_phases(gm, cpac) ≈ closure_phases(m, cpac)
    @test logclosure_amplitudes(gm, lcac) ≈ logclosure_amplitudes(m, lcac)


    rga = rand(gamp)
    rph = rand(gpha)
    @inferred logdensityof(gamp, rga)
    @inferred logdensityof(gpha, rph)
end
