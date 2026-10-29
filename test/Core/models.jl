using ChainRulesTestUtils
using ChainRulesCore
using FiniteDifferences
using Enzyme
using PythonCall
using FFTW
using StructArrays
using StaticArrays
using Distributions
using VLBIImagePriors
using StaticArrays
using Tables
import TransformVariables as TV
using VLBIImagePriors
using JET

ntequal(x::NamedTuple{N}, y::NamedTuple{N}) where {N} = map(_ntequal, (x), (y))
ntequal(x, y) = false
_ntequal(x::T, y::T) where {T <: NamedTuple} = ntequal(values(x), values(y))
_ntequal(x::T, y::T) where {T <: Tuple} = map(_ntequal, x, y)
_ntequal(x, y) = x ≈ y

function isequalmissing(x, y)
    xm = x |> ismissing |> collect
    ym = y |> ismissing |> collect
    return xm == ym
end


function test_caltable(c1, sites)
    @test Tables.istable(typeof(c1))
    @test Tables.rowaccess(typeof(c1))
    @test Tables.rows(c1) === c1
    @test Tables.columnaccess(c1)
    clmns = Tables.columns(c1)
    @test clmns[1] == Comrade.times(c1)
    @test clmns[2] == Comrade.frequencies(c1)
    @test Bool(prod(skipmissing(Tables.matrix(clmns)[:, (begin + 2):end]) .== skipmissing(Comrade.gmat(c1))))
    @test c1.Ti == Comrade.times(c1)
    @test c1.Ti == Tables.getcolumn(c1, 1)
    @test c1.Fr == Comrade.frequencies(c1)
    @test c1.Fr == Tables.getcolumn(c1, 2)
    @test isequalmissing(c1.AA, Tables.getcolumn(c1, 3))

    @test maximum(abs, skipmissing(c1.AA) .- skipmissing(Tables.getcolumn(c1, :AA))) ≈ 0
    @test maximum(abs, skipmissing(c1.AA) .- skipmissing(Tables.getcolumn(c1, 3))) ≈ 0
    @test Tables.columnnames(c1) == [:Ti, :Fr, sites...]

    c1row = Tables.getrow(c1, 30)
    @test eltype(c1) == typeof(c1row)
    @test c1row.Ti == c1.Ti[30]
    @test c1row.AA == c1.AA[30]
    @test c1row.Fr == c1.Fr[30]
    @test Tables.getcolumn(c1row, :AA) == c1.AA[30]
    @test Tables.getcolumn(c1row, :Ti) == c1.Ti[30]
    @test Tables.getcolumn(c1row, :Fr) == c1.Fr[30]
    @test Tables.getcolumn(c1row, 3) == c1.AA[30]
    @test Tables.getcolumn(c1row, 2) == c1.Fr[30]
    @test Tables.getcolumn(c1row, 1) == c1.Ti[30]
    @test propertynames(c1) == propertynames(c1row) == [:Ti, :Fr, sites...]
    @test Tables.getcolumn(c1row, Float64, 1, :Ti) == c1.Ti[30]
    @test Tables.getcolumn(c1row, Float64, 2, :Fr) == c1.Fr[30]
    @test Tables.getcolumn(c1row, Float64, 3, :AA) == c1.AA[30]
    @test isequalmissing(c1[1:10, :AA], c1.AA[1:10])
    @test isequalmissing(c1[[1, 2], :AA], c1.AA[[1, 2]])
    @test isequalmissing(@view(c1[1:10, :AA]), @view(c1.AA[1:10]))
    @test isequalmissing(@view(c1[[1, 2], :AA]), @view(c1.AA[[1, 2]]))

    Tables.schema(c1) isa Tables.Schema
    Tables.getcolumn(c1, Float64, 1, :test)
    Tables.getcolumn(c1, Float64, 2, :test)

    c1[1, :AA]
    c1[!, :AA]
    c1[:, :AA]
    @test length(c1) == length(c1.AA)
    @test c1[1, :] isa Comrade.CalTableRow
    @test length(Tables.getrow(c1, 1:5)) == 5

    @test plotcaltable(c1) isa Makie.Figure

    return show(c1)
end

@testset "SkyModel" begin

    f = test_model
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(f, test_prior(), g)
    show(IOBuffer(), MIME"text/plain"(), skym)

    x = rand(Comrade.NamedDist(test_prior()))
    m = Comrade.skymodel(skym, x)
    skyf = FixedSkyModel(m, g)

    @testset "ObservedSkyModel" begin
        _, vis, amp, lcamp, cphase = load_data()

        oskym, = Comrade.set_array(skym, vis)
        oskyf, = Comrade.set_array(skyf, vis)

        @test Comrade.skymodel(oskym, x) == m
        dmimg, dmvis = Comrade.idealmaps(Comrade.VisData(), oskym, (; sky = x))
        dfimg, dfvis = Comrade.idealmaps(Comrade.VisData(), oskyf, (; sky = x))
        @test dmimg ≈ dfimg
        @test dmvis ≈ dfvis
    end

    @testset "MultiSkyModel" begin
        _, vis, amp, lcamp, cphase = load_data()
        skytot = MultiSkyModel((dynamic = skym, static = skyf))
        oskytot, ptot = Comrade.set_array(skytot, vis)

        show(IOBuffer(), MIME"text/plain"(), skytot)


        x = rand(ptot)

        oskym, = Comrade.set_array(skym, vis)
        oskyf, = Comrade.set_array(skyf, vis)

        tt = Comrade.skymodel(oskytot, x)
        @test tt.dynamic == Comrade.skymodel(oskym, x.dynamic)
        @test tt.static == Comrade.skymodel(oskyf, x.static)

        vtot = last(Comrade.idealmaps(Comrade.VisData(), oskytot, (; sky = x)))
        vdyn = last(Comrade.idealmaps(Comrade.VisData(), oskym, (; sky = x.dynamic)))
        vstat = last(Comrade.idealmaps(Comrade.VisData(), oskyf, (; sky = x.static)))

        @test_opt Comrade.idealmaps(Comrade.VisData(), oskytot, (; sky = x))

        @test vtot ≈ vdyn + vstat
    end
end

@testset "GMRF" begin
    δ = rand(GMRF(64.0, (128, 128)))
    m = modify(Gaussian(), Stretch(μas2rad(40.0)))
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 128, 128)
    img = apply_fluctuations(exp, m, g, δ)
    mimg = intensitymap(m, g)
    img2 = apply_fluctuations(exp, mimg, δ)
    @test img isa IntensityMap
    @test img2 ≈ img

    img3 = apply_fluctuations(mimg, δ)

    img4 = apply_fluctuations(CenteredLR(), m, g, δ)
    img5nf = apply_fluctuations(CenteredLR(), mimg, δ)
    img5 = apply_fluctuations(CenteredLR(), mimg ./ Comrade.flux(mimg), δ)
    @test img4 ≈ img5
    @test img5nf ≈ img5

    img6 = apply_fluctuations(UnitFluxMap(exp), mimg ./ Comrade.flux(mimg), δ)
    @test img6 ≈ img5
    img7 = apply_fluctuations(UnitFluxMap(exp), mimg, δ)
    @test img7 ≈ img5
end

function FiniteDifferences.to_vec(k::SiteArray)
    v, b = to_vec(parent(k))
    back(x) = SiteArray(b(x), k.times, k.frequencies, k.sites)
    return v, back
end


@testset "InstrumentModel" begin
    _, dvis, amp, lcamp, cphase, dcoh = load_data()

    @testset "rotation" begin
        @test rotation(0.3, CirBasis()) ≈ SMatrix{2, 2}(cis(-0.3), 0, 0, cis(0.3))
        for B in (CirBasis(), LinBasis())
            R = rotation(0.3, B)
            @test R' * R ≈ I
            @test rotation(0.1, B) * rotation(0.2, B) ≈ R
            Tc = basis_transform(CirBasis(), B)
            @test R ≈ Tc * rotation(0.3, CirBasis()) * Tc'
        end
    end

    @testset "SiteArray" begin

        G = Jones(
            (x, c) -> exp(x.lg + 1im * x.gp);
            lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
        )

        intm = InstrumentModel(G)
        ointm, printm = bind_instrument(intm, dvis)
        show(IOBuffer(), MIME"text/plain"(), ointm)
        x = rand(printm)
        sl = Comrade.SiteLookup(x.lg)

        @test sl isa Comrade.SiteLookup
        s2 = SiteArray(parent(x.lg), sl)
        @test s2 == x.lg
        @test Comrade.times(s2) == Comrade.times(x.lg)
        @test Comrade.frequencies(s2) == Comrade.frequencies(x.lg)
        @test Comrade.sites(s2) == Comrade.sites(x.lg)


        @test x.lg[1] == parent(x.lg)[1]
        x.lg[1] = 1.0
        @test x.lg[1] == 1.0

        sarr = x.lg[1:16]
        @test_throws DimensionMismatch similar(sarr, Float64, (4, 4)) isa Comrade.SiteArray

        @test Comrade.SiteArray(x.lg, sl) == x.lg

        @inferred Comrade.time(x.lg, 5.0 .. 6.0)
        @inferred Comrade.frequency(x.lg, 1.0 .. 400.0)

        @test x.lg ≈ SiteArray(x.lg, x.lg.times, x.lg.frequencies, x.lg.sites)
        @inferred x.lg[1, 1, 1]
        x.lg[1, 1, 1, 1] = 1.0
        @test x.lg[1] ≈ 1.0

        # ps = ProjectTo(x.lg)
        # @test ps(x.lg) == x.lg
        # @test ps(NoTangent()) isa NoTangent
        # @test ps(Tangent{typeof(x.lg)}(data = parent(x.lg))) == x.lg

    end

    @testset "StokesI" begin
        t = dvis
        vis = parent(Comrade.measurement(t))

        G = Jones(
            (x, c) -> exp(x.lg + 1im * x.gp);
            lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
        )

        intm = InstrumentModel(G)
        show(IOBuffer(), MIME"text/plain"(), intm)

        ointm, printm = bind_instrument(intm, t)
        x = rand(printm)
        x.lg .= 0
        x.gp .= 0
        vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
        # test_rrule(Comrade.apply_instrument, vis, ointm⊢NoTangent(), (;instrument=x))
        @test vout ≈ vis


        ointid, pr = bind_instrument(Comrade.IdealInstrumentModel(), t)
        vout = Comrade.apply_instrument(vis, ointid, (;))
        @test vout ≈ vis

        # Now check that everything is being applied right
        bl = baselines(t)
        for s in sites(t)
            x.lg .= 0
            x.gp .= 0

            inds1 = findall(b -> (b[1] == s), bl)
            inds2 = findall(b -> (b[2] == s), bl)
            ninds = findall(b -> (b[1] != s && b[2] != s), bl)

            # Now amps
            x.lg .= 0
            xlgs = x.lg[S = s]
            xlgs .= log(2)
            vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
            @test vout[inds1, :] ≈ 2 .* vis[inds1, :]
            @test vout[inds2, :] ≈ 2 .* vis[inds2, :]
            @test vout[ninds, :] ≈ vis[ninds, :]

            # Now Phases
            x.lg .= 0
            xgps = x.gp[S = s]
            xgps .= π / 4
            vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
            @test vout[inds1, :] ≈ vis[inds1, :] .* exp(1im * π / 4)
            @test vout[inds2, :] ≈ vis[inds2, :] .* exp(-1im * π / 4)
            @test vout[ninds, :] ≈ vis[ninds, :]

            # Now Phases and amps
            x.lg .= 0
            x.gp .= 0
            xlgs = x.lg[S = s]
            xlgs .= log(2)
            xgps = x.gp[S = s]
            xgps .= π / 4
            vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
            @test vout[inds1, :] ≈ vis[inds1, :] .* exp(log(2) + 1im * π / 4)
            @test vout[inds2, :] ≈ vis[inds2, :] .* exp(log(2) - 1im * π / 4)
            @test vout[ninds, :] ≈ vis[ninds, :]
        end

    end


    @testset "Multifrequency parameters are shared by the channels" begin
        ν0 = only(info(dvis).channels.frequency)
        dvismf = with_channel(dvis, ν0 / 2)
        ts = Comrade.timestamps(IntegSeg(), dvismf)
        @test all(t -> t.dt > 0, ts)

        gp = ArrayPrior(
            IIDSitePrior(IntegSeg(), DiagonalVonMises(0.0, inv(π^2)));
            refant = SensitiveReference(0.0),
        )
        smap = Comrade.build_sitemap(gp, dvismf)
        @test length(unique(smap.frequencies)) == 1
        @test smap.sites == Comrade.build_sitemap(gp, dvis).sites
        @test Comrade.ObservedArrayPrior(gp, dvismf) isa Comrade.ObservedArrayPrior

        # A delay: the gain phase grows linearly with frequency from one parameter per site.
        D = Jones(
            (x, c) -> exp(2π * 1im * x.τ * (c.Fr - ν0));
            τ = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0e-10))),
        )
        ointm, printm = bind_instrument(InstrumentModel(D), dvismf)
        x = rand(printm)
        x.τ .= 0
        s = first(sites(dvismf))
        xτs = x.τ[S = s]
        xτs .= 1.0e-10
        vis = parent(Comrade.measurement(dvismf))
        vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
        bl = baselines(dvismf)
        i1 = findall(b -> b[1] == s, bl)
        ν = info(dvismf).channels.frequency
        @test length(x.τ) == length(Comrade.build_sitemap(D.priors.τ, dvis).sites)
        @test vout[i1, :] ≈ vis[i1, :] .* exp.(2π * 1im * 1.0e-10 .* (ν' .- ν0))
    end


    @testset "Coherencies" begin
        tcoh = dcoh
        vis = StokesParams.(CoherencyMatrix.(collect(parent(Comrade.measurement(tcoh))), Ref(CirBasis())))
        G = Jones(;
            lgR = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            gpR = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
            lgrat = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            gprat = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)), refant = SingleReference(:AA, 0.0)),
        ) do x, c
            gR = exp(x.lgR + 1im * x.gpR)
            gL = gR * exp(x.lgrat + 1im * x.gprat)
            return gains(gR, gL)
        end

        D = Jones(;
            dRx = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2))),
            dRy = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2))),
            dLx = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2))),
            dLy = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2))),
        ) do x, c
            dR = complex(x.dRx, x.dRy)
            dL = complex(x.dLx, x.dLy)
            return leakage(dR, dL)
        end

        J = G * D
        # A nested product flattens to the same factors.
        J2 = G * (D * Comrade.JonesProduct(()))

        JG = Jones((x, c) -> SMatrix{2, 2}(x.lg, x.lg, x.lg, x.lg); lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))))


        intm = InstrumentModel(J)
        intm2 = InstrumentModel(J2)
        intjg = InstrumentModel(JG)
        show(IOBuffer(), MIME"text/plain"(), intm)


        ointm, printm = bind_instrument(intm, tcoh)
        ointm2, printm2 = bind_instrument(intm2, tcoh)
        ointjg, printjg = bind_instrument(intjg, tcoh)

        x = rand(printjg)
        fj = Comrade.forward_jones(ointjg, (; instrument = x), tcoh)
        L = Comrade.layout(ointjg)
        U = [Comrade._load(L.response, k) for k in eachindex(L.coords)]
        c = first(L.coords)
        l = parent(x.lg)[L.entries.lg[1]]
        @test fj[c.ant, findfirst(==(c.Ti), Comrade.times(tcoh)), 1] ≈ SMatrix{2, 2}(l, l, l, l) * U[1]

        # In the circular basis the response is the feed rotation alone.
        st = info(tcoh).stations
        φ = [feed_rotation(st.mount[c.ant], c.parallactic, c.elevation, c.azimuth) + st.feed_angles[c.ant][1] for c in L.coords]
        @test collect(U) ≈ map(a -> SMatrix{2, 2}(cis(-a), 0, 0, cis(a)), φ)

        @testset "ObservedArrayPrior" begin
            @inferred logpdf(printm, rand(printm))
            @inferred logpdf(printm2, rand(printm2))
            x = rand(printm)
            @test logpdf(printm, x) ≈ logpdf(printm2, x)
            @test Comrade.transport_node(printm, Comrade.TVFlat()) isa TV.AbstractTransform
            p = rand(printm)
            t = Comrade.transport_node(printm, Comrade.TVFlat())
            pout = TV.transform(t, TV.inverse(t, p))
            dp = ntequal(p, pout)
            @test dp.lgR
            @test dp.lgrat
            @test dp.gprat
            @test dp.dRx
            @test dp.dRy
            @test dp.dLx
            @test dp.dLy
        end

        @testset "ObservedArrayPrior std-space transports" begin
            # The StdUniform/StdNormal (cube/normal) path uses the ProbabilityTransports
            # node interface (`pfwd_step`/`pback_step!`, Jacobian-free), as opposed to the
            # TVFlat TransformVariables path tested above. Verify the forward/backward
            # transports round-trip a parameter draw.
            for space in (Comrade.StdUniform(), Comrade.PT.StdNormal())
                t = Comrade.transport_node(printm, space)
                p = rand(printm)
                u = Comrade.latent_pback(t, p)
                @test length(u) == Comrade.PT.dimension(t)
                pout = Comrade.latent_pfwd(t, u)
                dp = ntequal(p, pout)
                @test dp.lgR
                @test dp.gpR
                @test dp.lgrat
                @test dp.gprat
                @test dp.dRx
                @test dp.dRy
                @test dp.dLx
                @test dp.dLy
            end
        end
        pintm, _ = bind_instrument(InstrumentModel(), tcoh)


        x = rand(printm)
        x.lgR .= 0
        x.lgrat .= 0
        x.gpR .= 0
        x.gprat .= 0
        x.dRx .= 0
        x.dRy .= 0
        x.dLx .= 0
        x.dLy .= 0

        vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
        vper = Comrade.apply_instrument(vis, pintm, (; instrument = NamedTuple()))
        @test vout ≈ vper

        # test_rrule(Comrade.apply_instrument, vis, ointm⊢NoTangent(), (;instrument=x))

        # # Now check that everything is being applied right
        bl = baselines(tcoh)
        for s in sites(tcoh)
            x.lgR .= 0
            x.lgrat .= 0
            x.gpR .= 0
            x.gprat .= 0
            x.dRx .= 0
            x.dRy .= 0
            x.dLx .= 0
            x.dLy .= 0


            inds1 = findall(b -> (b[1] == s), bl)
            inds2 = findall(b -> (b[2] == s), bl)
            ninds = findall(b -> (b[1] != s && b[2] != s), bl)

            # Now amp-offsets
            x.lgR .= 0
            x.lgrat .= 0
            x.gpR .= 0
            x.gprat .= 0
            x.dRx .= 0
            x.dRy .= 0
            x.dLx .= 0
            x.dLy .= 0

            xlgRs = x.lgR[S = s]
            xlgRs .= log(2)
            xlgrat = x.lgrat[S = s]
            xlgrat .= -log(2)
            vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
            G = SMatrix{2, 2}(2.0, 0.0, 0.0, 1.0)
            @test vout[inds1, :] ≈ Ref(G) .* vper[inds1, :]
            @test vout[inds2, :] ≈ vper[inds2, :] .* Ref(G)
            @test vout[ninds, :] ≈ vper[ninds, :]

            # Now phases
            x.lgR .= 0
            x.lgrat .= 0
            x.gpR .= 0
            x.gprat .= 0
            x.dRx .= 0
            x.dRy .= 0
            x.dLx .= 0
            x.dLy .= 0

            xgpRs = x.gpR[S = s]
            xgpRs .= π / 3
            xgprat = x.gprat[S = s]
            xgprat .= -π / 3
            vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
            G = SMatrix{2, 2}(exp(1im * π / 3), 0.0, 0.0, exp(1im * 0.0))
            @test vout[inds1, :] ≈ Ref(G) .* vper[inds1, :]
            @test vout[inds2, :] ≈ vper[inds2, :] .* Ref(adjoint(G))
            @test vout[ninds, :] ≈ vper[ninds, :]


            # Now dterms
            x.lgR .= 0
            x.lgrat .= 0
            x.gpR .= 0
            x.gprat .= 0
            x.dRx .= 0
            x.dRy .= 0
            x.dLx .= 0
            x.dLy .= 0

            xdRxs = x.dRx[S = s]
            xdRxs .= 0.1
            xdRys = x.dRy[S = s]
            xdRys .= 0.2
            xdLxs = x.dLx[S = s]
            xdLxs .= 0.3
            xdLys = x.dLy[S = s]
            xdLys .= 0.4

            vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
            D = SMatrix{2, 2}(1.0, 0.3 + 0.4im, 0.1 + 0.2im, 1.0)
            @test vout[inds1, :] ≈ Ref(D) .* vper[inds1, :]
            @test vout[inds2, :] ≈ vper[inds2, :] .* Ref(adjoint(D))
            @test vout[ninds, :] ≈ vper[ninds, :]
        end

        @testset "caltable test" begin
            c1 = caltable(x.lgR)
            test_caltable(c1, sort(sites(amp)))
        end

    end


    @testset "PerStation" begin
        tcoh = dcoh
        vis = StokesParams.(CoherencyMatrix.(collect(parent(Comrade.measurement(tcoh))), Ref(CirBasis())))
        g0 = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))
        prior = (;
            lgR = ArrayPrior(g0), gpR = ArrayPrior(g0; refant = SingleReference(:AA, 0.0)),
            lgL = ArrayPrior(g0), gpL = ArrayPrior(g0; refant = SingleReference(:AA, 0.0)),
        )
        δ = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.3)))
        qwp = (x, c) -> SMatrix{2, 2}(1.0 + 0im, 0, 0, cis(x.δ))
        G = Jones((x, c) -> gains(exp(complex(x.lgR, x.gpR)), exp(complex(x.lgL, x.gpL))); prior...)
        Q = Jones(qwp; δ)

        oplain, pplain = bind_instrument(InstrumentModel(G), tcoh)
        mper = InstrumentModel(PerStation(G; AA = G * Q))
        @test mper.prior == merge(prior, (; δ))
        oper, pper = bind_instrument(mper, tcoh)
        x = rand(pper)
        @test unique(sites(x.δ)) == [:AA]
        @test length(x.δ) == 1
        xg = Base.structdiff(x, NamedTuple{(:δ,)})
        @test length(rand(pplain).lgR) == length(x.lgR)

        vplain = Comrade.apply_instrument(vis, oplain, (; instrument = xg))
        @test Comrade.apply_instrument(vis, oper, (; instrument = merge(x, (; δ = zero(x.δ))))) ≈ vplain

        δAA = only(parent(x.δ))
        aa = findfirst(==(:AA), info(tcoh).stations.name)
        Qref = Jones((x, c) -> SMatrix{2, 2}(1.0 + 0im, 0, 0, c.ant == aa ? cis(δAA) : one(δAA)))
        oref, _ = bind_instrument(InstrumentModel(G * Qref), tcoh)
        vper = Comrade.apply_instrument(vis, oper, (; instrument = x))
        @test vper ≈ Comrade.apply_instrument(vis, oref, (; instrument = xg))
        @test !(vper ≈ vplain)

        fj = Comrade.forward_jones(oper, (; instrument = x), tcoh)
        @test size(fj, 1) == length(info(tcoh).stations.name)

        ti = dvis
        Gi = Jones((x, c) -> exp(complex(x.lgR, x.gpR)); prior[(:lgR, :gpR)]...)
        @test_throws "different types" bind_instrument(InstrumentModel(PerStation(Gi; AA = Gi * Q)), ti)
        @test_throws "not in the table" bind_instrument(InstrumentModel(PerStation(G; XX = G * Q)), tcoh)
        orphan = InstrumentModel(PerStation(G), merge(prior, (; δ)), CirBasis(), ReportGauge(), false, false)
        @test_throws "no station's term reads" bind_instrument(orphan, tcoh)
        δAAonly = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.3)); stations = (:AA,))
        @test_throws "δ has no entry for station" bind_instrument(InstrumentModel(G * Jones(qwp; δ = δAAonly)), tcoh)
        lgR2 = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.2)))
        G2 = Jones((x, c) -> exp(x.lgR); lgR = lgR2)
        @test_throws "different priors" InstrumentModel(G * G2)
        @test_throws "different priors" InstrumentModel(PerStation(G; AA = G2))
        @test_throws "whole instrument" PerStation(G) * Q
        @test_throws "cannot be nested" PerStation(G; AA = PerStation(G))
    end

    @testset "Coherencies Multifrequency" begin
        tsi = dcoh
        tmf = with_channel(dcoh, 345.0e9)
        @test size(tmf) == (size(tsi, 1), 2)
        vissi = StokesParams.(CoherencyMatrix.(collect(parent(Comrade.measurement(tsi))), Ref(CirBasis())))
        vismf = StokesParams.(CoherencyMatrix.(collect(parent(Comrade.measurement(tmf))), Ref(CirBasis())))
        G = Jones(;
            lgR = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            gpR = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
            lgrat = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            gprat = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
        ) do x, c
            gR = exp(x.lgR + 1im * x.gpR)
            gL = gR * exp(x.lgrat + 1im * x.gprat)
            return gains(gR, gL)
        end

        D = Jones(;
            dRx = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2))),
            dRy = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2))),
            dLx = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2))),
            dLy = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2))),
        ) do x, c
            dR = complex(x.dRx, x.dRy)
            dL = complex(x.dLx, x.dLy)
            return leakage(dR, dL)
        end

        J = G * D


        intm = InstrumentModel(J)
        show(IOBuffer(), MIME"text/plain"(), intm)

        ointsi, printsi = bind_instrument(intm, tsi)
        ointmf, printmf = bind_instrument(intm, tmf)
        @testset "Site lookup" begin
            trsi = Comrade.transport_node(printsi, Comrade.TVFlat())
            trmf = Comrade.transport_node(printmf, Comrade.TVFlat())
            # every channel reads the same parameters
            @test trmf.lgR.site_map.lookup == trsi.lgR.site_map.lookup
            l = length(trsi.lgR.site_map.frequencies)
            @test reduce(vcat, trmf.lgR.site_map.lookup) |> sort == 1:l
        end

        Lsi, Lmf = Comrade.layout(ointsi), Comrade.layout(ointmf)
        @test length(Lmf.coords) == 2 * length(Lsi.coords)
        # The single-channel cell with the station and time of each two-channel cell
        cellkey(L, c) = (L.coords.ant[c], L.coords.Ti[c])
        sicell = Dict(cellkey(Lsi, c) => c for c in eachindex(Lsi.coords))
        same = [sicell[cellkey(Lmf, c)] for c in eachindex(Lmf.coords)]
        for ch in 1:2
            @test same[Lmf.cell1[:, ch]] == Lsi.cell1[:, 1]
            @test same[Lmf.cell2[:, ch]] == Lsi.cell2[:, 1]
        end
        # The response does not depend on the channel.
        responses(L) = [Comrade._load(L.response, k) for k in eachindex(L.coords)]
        @test responses(Lmf) ≈ responses(Lsi)[same]
        # Each channel's cells read the entries the single-channel cells read.
        for ν in unique(Lmf.coords.Fr)
            chν = findall(==(ν), Lmf.coords.Fr)
            for p in propertynames(Lsi.entries)
                @test Lmf.entries[p][chν] == Lsi.entries[p][same[chν]]
            end
        end

        xsi = rand(printsi)
        xmf = rand(printmf)
        @test map(length, xmf) == map(length, xsi)
        foreach((a, b) -> copyto!(parent(a), parent(b)), xmf, xsi)
        vsi = Comrade.apply_instrument(vissi, ointsi, (; instrument = xsi))
        vmf = Comrade.apply_instrument(vismf, ointmf, (; instrument = xmf))
        @test vmf[:, 1] ≈ vsi[:, 1]
        @test vmf[:, 2] ≈ vsi[:, 1]

        @testset "caltable test" begin
            xmf = rand(printmf)
            c1 = caltable(xmf.lgR)
            test_caltable(c1, sort(sites(amp)))
        end

    end

    @testset "Hypercube Instrument" begin
        _, dvis, amp, lcamp, cphase, dcoh = load_data()
        gmap = (x, c) -> exp(x.lg + 1im .* x.gp)

        @testset "standard" begin
            intprior = (
                lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
                gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
            )

            intm = InstrumentModel(Jones(gmap; intprior...))
            skym = SkyModel(test_model, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 256, 256))
            post = VLBIPosterior(skym, intm, dvis)

            cpost = ascube(post)
            x = prior_sample(cpost)
            p = Comrade.transform(cpost, x)
            @test all(x -> 0 < x < 1, x)

            l1 = logdensityof(cpost, x)
            l2 = Comrade.loglikelihood(post, p)

            @test l1 ≈ l2
        end


    end

    @testset "Integration" begin
        _, dvis, amp, lcamp, cphase, dcoh = load_data()
        ts = Comrade.timestamps(ScanSeg(), dvis)
        tt = Comrade.timestamps(TrackSeg(), dvis)
        ti = Comrade.timestamps(IntegSeg(), dvis)
        @test length(tt) < length(ts) ≤ length(ti)
    end

    @testset "EpochSeg" begin
        _, dvis, = load_data()
        two = vcat(dvis, later(dvis, 1))
        tr = info(two).tracks
        @test length(tr) == 2
        @test length(Comrade.timestamps(EpochSeg(), dvis)) == 1
        @test length(Comrade.timestamps(TrackSeg(), two)) == 1
        te = Comrade.timestamps(EpochSeg(), two)
        @test length(te) == 2
        @test all(τ -> count(s -> τ ∈ s, te) == 1, domain(two).Ti)
        @test [s.t0 for s in te] ≈ (tr.start .+ tr.stop) ./ 2

        G = Jones((x, c) -> exp(x.lg); lg = ArrayPrior(IIDSitePrior(EpochSeg(), Normal(0.0, 0.1))))
        ointm, printm = bind_instrument(InstrumentModel(G), two)
        x = rand(printm)
        @test length(x.lg) == 2 * length(sites(dvis))
        parent(x.lg) .= map(s -> tr.start[2] ∈ s ? 0.1 : 0.0, Comrade.times(x.lg))
        vis = parent(measurement(two))
        vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))
        day2 = domain(two).Ti .> tr.start[2]
        @test vout[.!day2, :] ≈ vis[.!day2, :]
        @test vout[day2, :] ≈ vis[day2, :] .* exp(0.2)
    end

    @testset "IntegrationTime" begin
        ti = Comrade.IntegrationTime(10, 5.0, 0.1)
        @test Comrade.mjd(ti) == ti.mjd
        @test ti.t0 ∈ Comrade.interval(ti)
        @test Comrade._center(ti) == ti.t0
        @test Comrade._region(ti) == 0.1
    end

    @testset "FrequencyChannel" begin
        fc = Comrade.FrequencyChannel(230.0e9, 8.0e9, 1)
        @test Comrade.channel(fc) == 1
        @test fc.central ∈ Comrade.interval(fc)
        @test Comrade._center(fc) == fc.central
        @test Comrade._region(fc) == 8.0e9
        @test 86.0e9 < fc
        @test fc < 345.0e9
    end


end
