using Distributions
using Plots
using Tables

@testset "gains" begin
    _,vis, amp, lcamp, cphase, = load_data()

    st = scantable(vis)
    stamp = scantable(amp)
    stcp = scantable(cphase)
    stlca = scantable(lcamp)


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
    gamp_prior = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))
    gph_prior_0 = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.5), length(tel)))
    gph_prior = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))

    jcache = jonescache(vis, ScanSeg())
    jcacher = jonescache(vis, ScanSeg(), true)

    # test the design matrix
    d1 = Comrade.DesignMatrix(jcache.m1, jcache.times, jcache.stations)
    @test size(d1) == size(jcache.m1)
    @test d1[1,1] == jcache.m1[1,1]
    @test Base.IndexStyle(d1) == Base.IndexStyle(jcache.m1)

    similar(d1)


    gamp = CalPrior(gamp_prior, jcache)
    gpha = CalPrior(gph_prior,  jcache)
    gphar = CalPrior(gph_prior_0, gph_prior, jcacher)

    ga = fill(1.0, size(rand(gamp))...)
    gm = JonesModel(jonesStokes(ga, jcache), m)

    @inferred Comrade.intensity_point(gm, (X=0.0, Y=0.0))
    @test gm.model === m

    g = imagepixels(μas2rad(200.0), μas2rad(200.0), 128, 128)
    img1 = intensitymap(gm, g)
    img2 = intensitymap(m, g)

    @test img1 ≈ img2

    @test flux(m) ≈ flux(gm)

    j1 = jonesStokes(exp.(ga), jcache)
    j2 = jonesStokes(exp, ga, jcache)

    @test j1.m1 == j2.m1
    @test j1.m2 == j2.m2

    c1 = caltable(jcache, ga)
    c1r = caltable(jcacher, ga)
    @test c1.time == c1r.time
    @test cumsum(c1.AA) == c1r.AA
    c2 = caltable(amp, ga)

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


    @test prod(skipmissing(Comrade.gmat(c2) .== Comrade.gmat(c1)))

    plot(c1, layout=(3,3), size=(600,500))

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

@testset "dterms" begin
    _,vis, amp, lcamp, cphase, dcoh = load_data()
    tel = stations(vis)

    dReal = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))
    dImag = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))


    dcache = jonescache(dcoh, TrackSeg())
    pdReal = CalPrior(dReal, dcache)
    pdImag = CalPrior(dImag, dcache)
    dRx = rand(pdReal)
    dRy = rand(pdReal)
    jd1 = jonesD(identity, dRx, dRy, dcache)
    jd2 = jonesD(dRx,  dRy, dcache)

    @test jd1.m1 == jd2.m1
    @test jd1.m2 == jd2.m2

    @inferred jonesD(exp, rand(pdReal), rand(pdImag), dcache)

    caltable(dcache, complex.(dRx, dRy))

end


@testset "rlgains" begin
    _,vis, amp, lcamp, cphase, dcoh = load_data()
    tel = stations(vis)

    gprior0 = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.5), length(tel)))
    gprior1 = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))

    trackcache = jonescache(dcoh, TrackSeg())
    scancache = jonescache(dcoh, ScanSeg())

    dlgp = CalPrior(gprior0, scancache)
    dgpp = CalPrior(gprior1, scancache)
    dlgr = CalPrior(gprior0, scancache)
    dgpr = CalPrior(gprior1,scancache)

    lgp = rand(dlgp)
    gpp = rand(dgpp)
    lgr = rand(dlgr)
    gpr = rand(dgpr)

    gpro = exp.(lgp .+ 1im.*gpp)
    grat = exp.(lgr .+ 1im.*gpr)
    Gp = jonesG(gpro, gpro, scancache)
    Gr = jonesG(grat, inv.(grat), scancache)

    j1 = Gp*Gr
    j2 = jonesG(gpro.*grat, gpro./grat, scancache)

    @test j1.m1 ≈ j2.m1
    @test j1.m2 ≈ j2.m2

    m1 = map(x->getproperty(x, :m1), (Gp,Gr))
    m2 = map(x->getproperty(x, :m2), (Gp,Gr))

    test_rrule(Comrade._allmul, m1, m2)

end
