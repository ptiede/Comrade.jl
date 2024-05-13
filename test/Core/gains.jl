using Distributions
using Plots
using Tables
using StructArrays
using StaticArrays

@testset "gains" begin
    _,vis, amp, lcamp, cphase, = load_data()

    st = timetable(vis)
    stamp = timetable(amp)
    stcp = timetable(cphase)
    stlca = timetable(lcamp)


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

    tel = sites(vis)
    gamp_prior = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))
    gph_prior_0 = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.5), length(tel)))
    gph_prior = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))

    jcache = jonescache(vis, ScanSeg())
    jcacher = jonescache(vis, ScanSeg(true))

    @test site_tuple(tel, ScanSeg()) == site_tuple(vis, ScanSeg())
    @test site_tuple(tel, ScanSeg(); AA=FixedSeg(0.0)) == site_tuple(vis, ScanSeg(); AA=FixedSeg(0.0))

    # test the design matrix
    d1 = Comrade.DesignMatrix(jcache.m1, jcache.schema.times, jcache.schema.sites)
    v = rand(size(d1, 2))
    @test d1*v ≈ jcache.m1*v
    @test size(d1) == size(jcache.m1)
    @test d1[1,1] == jcache.m1[1,1]
    @test Base.IndexStyle(d1) == Base.IndexStyle(jcache.m1)

    similar(d1)


    gamp = CalPrior(gamp_prior, jcache)
    gpha = CalPrior(gph_prior,  jcache)
    gphar = CalPrior(gph_prior_0, gph_prior, jcache)


    # gamp_h = HierarchicalCalPrior{Normal}(gamp_prior, site_tuple(vis, truncated(Normal(0.0, 0.1); lower=0.0)), jcache)
    # x = rand(gamp_h)
    # @inferred logdensityof(gamp_h, x)

    ga = fill(1.0, size(rand(gamp))...)
    gm = Comrade.VLBIModel(JonesModel(jonesStokes(ga, jcache)), m)

    @inferred Comrade.intensity_point(gm, (X=0.0, Y=0.0))
    @test gm.sky === m

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
        @test Tables.columnnames(c1) == [:time, sort(sites(amp))...]

        c1row = first(c1)
        @test eltype(c1) == typeof(c1row)
        @test c1row.time == c1.time[1]
        @test c1row.AA == c1.AA[1]
        @test Tables.getcolumn(c1row, :AA) == c1.AA[1]
        @test Tables.getcolumn(c1row, :time) == c1.time[1]
        @test Tables.getcolumn(c1row, 2) == c1.AA[1]
        @test Tables.getcolumn(c1row, 1) == c1.time[1]
        @test propertynames(c1) == propertynames(c1row) == [:time, sort(sites(amp))...]


    end


    @test prod(skipmissing(Comrade.gmat(c2) .== Comrade.gmat(c1)))

    plot(c1, layout=(3,3), size=(600,500))

    ac = arrayconfig(vis)
    cpac = arrayconfig(cphase)
    lcac = arrayconfig(lcamp)
    @test visibilitymap(gm, ac) ≈ visibilitymap(m, ac)
    @test amplitudes(gm, ac) ≈ amplitudes(m, ac)
    @test closure_phases(gm, cpac) ≈ closure_phases(m, cpac)
    @test logclosure_amplitudes(gm, lcac) ≈ logclosure_amplitudes(m, lcac)


    rga = rand(gamp)
    rph = rand(gpha)
    @inferred logdensityof(gamp, rga)
    @inferred logdensityof(gpha, rph)
end

@testset "calibration priors" begin
    _,vis, amp, lcamp, cphase, dcoh = load_data()
    tel = sites(vis)

    dReal = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))
    dImag = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))


    dcache = jonescache(dcoh, ScanSeg())
    dcache_si = jonescache(dcoh, ScanSeg(); autoref=SingleReference(:AA, 1.0+0.0im))
    pdReal = CalPrior(dReal, dcache)
    dcache_ra = jonescache(dcoh, ScanSeg(); autoref=RandomReference(1.0+0.0im))
    pdRA = CalPrior(dReal, dcache_ra)
    dcache_se0 = jonescache(dcoh, ScanSeg(); autoref=SEFDReference(1.0+0.0im, 0))
    dcache_se1 = jonescache(dcoh, ScanSeg(); autoref=SEFDReference(1.0+0.0im, 1))
    pdSE_0 = CalPrior(dReal, dcache_se0)
    pdSE_1 = CalPrior(dReal, dcache_se1)

    @test mean(pdReal) ≈ zeros(length(pdReal))
    @test var(pdReal) ≈ fill(0.1^2, length(pdReal))


    asflat(pdReal)
    asflat(pdRA)
    ascube(pdSE_0)
    ascube(pdSE_1)
end

@testset "dterms" begin
    _,vis, amp, lcamp, cphase, dcoh = load_data()
    tel = sites(vis)

    dReal = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))
    dImag = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))

    dReal = site_tuple(dcoh, Uniform(0.0, 1.0), AA = Uniform(-1.0, 0.0))
    dImag = site_tuple(dcoh, Uniform(0.0, 1.0), AA = Uniform(-1.0, 0.0))

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
    n = length(dcoh)

    visS = StructArray{StokesParams{Float64}}((
        fill(complex(0.0), n),
        fill(complex(2.0), n),
        fill(complex(0.0), n),
        fill(complex(0.0), n),
            ))

    vis = CoherencyMatrix.(visS, Ref(CirBasis()))

    function f(Q, d1, d2)
        dR1 = d1[1,2]
        dR2 = d2[1,2]
        dL1 = d1[2,1]
        dL2 = d2[2,1]
        return Q*SMatrix{2,2}(dR2' + dR1, dL1*dR2' + 1, 1+dR1*dL2', dL1 + dL2')
    end


    dR_1 = [if d > 0; 1.0; else -1.0 end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    dL_1 = [if d > 0; 1.0; else -1.0 end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    jd = jonesD(dR_1, dL_1, dcache)
    c = Comrade.apply_instrument(visS, JonesModel(jd))
    @test c == f.(2.0, jd.m1, jd.m2)

    dR_1 = [if d > 0; 1.0im; else -1.0im end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    dL_1 = [if d > 0; 1.0im; else 1.0im end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    jd = jonesD(dR_1, dL_1, dcache)
    c = Comrade.apply_instrument(visS, JonesModel(jd))
    @test c == f.(2.0, jd.m1, jd.m2)


    dR_1 = [if d > 0; 1.0im; else -1.0im end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    dL_1 = [if d > 0; 1.0+0.0im; else 1.0im end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    jd = jonesD(dR_1, dL_1, dcache)
    c = Comrade.apply_instrument(visS, JonesModel(jd))
    @test c == f.(2.0, jd.m1, jd.m2)


    visS = StructArray{StokesParams{Float64}}((
        fill(complex(0.0), n),
        fill(complex(0.0), n),
        fill(complex(0.0), n),
        fill(complex(2.0), n),
            ))

    vis = CoherencyMatrix.(visS, Ref(CirBasis()))

    function f(V, d1, d2)
        dR1 = d1[1,2]
        dR2 = d2[1,2]
        dL1 = d1[2,1]
        dL2 = d2[2,1]
        return V*SMatrix{2,2}(1 - dR1*dR2', dL1-dR2', dL2' - dR1, dL1*dL2' - 1)
    end


    dR_1 = [if d > 0; 1.0; else -1.0 end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    dL_1 = [if d > 0; 1.0; else -1.0 end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    jd = jonesD(dR_1, dL_1, dcache)
    c = Comrade.apply_instrument(visS, JonesModel(jd))
    @test c == f.(2.0, jd.m1, jd.m2)

    dR_1 = [if d > 0; 1.0im; else -1.0im end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    dL_1 = [if d > 0; 1.0im; else 1.0im end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    jd = jonesD(dR_1, dL_1, dcache)
    c = Comrade.apply_instrument(visS, JonesModel(jd))
    @test c == f.(2.0, jd.m1, jd.m2)


    dR_1 = [if d > 0; 1.0im; else -1.0im end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    dL_1 = [if d > 0; 1.0+0.0im; else 1.0im end for d in dRx]# dRx.*0.0 .+ 1.0 .+ 0.0im
    jd = jonesD(dR_1, dL_1, dcache)
    c = Comrade.apply_instrument(visS, JonesModel(jd))
    @test c == f.(2.0, jd.m1, jd.m2)

    test_rrule(Comrade._apply_instrument, vis, jd.m1, jd.m2, CirBasis()⊢NoTangent())

end
