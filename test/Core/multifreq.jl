using Distributions
using VLBIImagePriors
using Test

# Multifrequency instrument-prior structure tests. `build_mfvis` is defined in
# models.jl, which core.jl includes before this file.
@testset "Multifrequency site chains" begin
    _, dvis, amp, lcamp, cphase, dcoh = load_data()

    dvis2 = deepcopy(dvis)
    dvis2.config[:Fr] .= 345.0e9
    dvismf = build_mfvis(dvis, dvis2)
    arr1 = arrayconfig(dvis)
    arr2 = arrayconfig(dvis2)
    arrmf = arrayconfig(dvismf)

    # A ragged variant: drop :AA entirely from the second channel (GMVA-like coverage)
    keep = findall(bl -> !(:AA in bl), Comrade.datatable(arrayconfig(dvis2)).sites)
    dvis2r = Comrade.EHTObservationTable{Comrade.datumtype(dvis2)}(
        Comrade.measurement(dvis2)[keep], Comrade.noise(dvis2)[keep],
        Comrade.EHTArrayConfiguration(
            dvis2.config.bandwidth, dvis2.config.tarr, dvis2.config.scans,
            dvis2.config.mjd, dvis2.config.ra, dvis2.config.dec, dvis2.config.source,
            :UTC, Comrade.StructArray(Comrade.datatable(arrayconfig(dvis2))[keep])
        )
    )
    dvismfr = build_mfvis(dvis, dvis2r)
    arrmfr = arrayconfig(dvismfr)

    # independent brute-force point matcher against a lookup's segment edges (reference
    # implementation for the baselinemap and param-lookup tests)
    inseg(sl, i, t, f, s) =
        (sl.tlo[sl.tcode[i]] ≤ t < sl.thi[sl.tcode[i]]) &&
        (sl.flo[sl.fcode[i]] ≤ f < sl.fhi[sl.fcode[i]]) &&
        (sl.saxis[sl.scode[i]] == s)

    @testset "chain structure" begin
        pr = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        smap = Comrade.ObservedArrayPrior(pr, arrmf).sitemap
        smapr = Comrade.ObservedArrayPrior(pr, arrmfr).sitemap

        cs = Comrade.sitechains(smap)
        # every site observed both channels, ascending in frequency
        @test all(chs -> length(chs) == 2, values(cs))
        @test all(chs -> chs[1].fcode < chs[2].fcode, values(cs))
        @test all(chs -> Comrade.chainfreq(smap, chs[1]) < Comrade.chainfreq(smap, chs[2]), values(cs))
        # chain indices are ascending in time and cover the storage exactly once
        @test all(chs -> all(c -> issorted(c.tcodes), chs), values(cs))
        allinds = sort(reduce(vcat, [c.inds for chs in values(cs) for c in chs]))
        @test allinds == 1:length(smap.tcode)

        # ragged coverage: :AA only has the low channel, everyone else has both
        csr = Comrade.sitechains(smapr)
        @test length(csr.AA) == 1
        @test csr.AA[1].fcode == 1
        @test all(s -> s == :AA || length(csr[s]) == 2, keys(csr))
        allindsr = sort(reduce(vcat, [c.inds for chs in values(csr) for c in chs]))
        @test allindsr == 1:length(smapr.tcode)
    end

    @testset "multifreq GM logpdf = sum of per-channel chains" begin
        for anchored in (false, true)
            pr = ArrayPrior(
                GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0); anchored);
                refant = SEFDReference(0.0)
            )
            omf = Comrade.ObservedArrayPrior(pr, arrmf)
            o1 = Comrade.ObservedArrayPrior(pr, arr1)
            o2 = Comrade.ObservedArrayPrior(pr, arr2)

            # mirror configs: storage order is freq-outer, so the multifreq vector is the
            # concatenation of the per-channel vectors
            @test length(omf) == length(o1) + length(o2)
            x1 = rand(o1)
            x2 = rand(o2)
            xmf = vcat(parent(x1), parent(x2))
            @test logpdf(omf.dists, xmf) ≈ logpdf(o1.dists, parent(x1)) + logpdf(o2.dists, parent(x2))
        end
    end

    @testset "hierarchical multifreq shares site hyperparameters" begin
        pr = ArrayPrior(
            GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(1.0), τ = Exponential(2.0)));
            refant = SEFDReference(0.0)
        )
        omf = Comrade.ObservedArrayPrior(pr, arrmf)
        o1 = Comrade.ObservedArrayPrior(pr, arr1)
        o2 = Comrade.ObservedArrayPrior(pr, arr2)
        s = rand(omf)
        # one hyperparameter set per real site, no pseudo-sites
        @test sort(collect(keys(s.hyperparams))) ==
            sort(collect(keys(Comrade.sitechains(omf.sitemap))))
        # chain density factorizes over channels at shared hyperparameters
        hp = s.hyperparams
        x1 = rand(o1)
        x2 = rand(o2)
        xmf = vcat(parent(x1.params), parent(x2.params))
        @test Comrade._chain_logpdf(omf.dists, xmf, hp) ≈
            Comrade._chain_logpdf(o1.dists, parent(x1.params), hp) +
            Comrade._chain_logpdf(o2.dists, parent(x2.params), hp)
    end

    @testset "ragged transforms round-trip" begin
        for pr in (
                ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)); refant = SEFDReference(0.0)),
                ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0));
                    refant = SEFDReference(0.0)
                ),
            )
            obs = Comrade.ObservedArrayPrior(pr, arrmfr)
            t = Comrade.transport_node(obs, Comrade.TVFlat())
            n = Comrade.TV.dimension(t)
            x0 = 0.3 .* sin.(1:n)
            y = Comrade.TV.transform(t, x0)
            xinv = similar(x0)
            Comrade.TV.inverse_at!(xinv, firstindex(xinv), t, y)
            @test xinv ≈ x0
            @test isfinite(logpdf(obs.dists, parent(y)))
        end
    end

    @testset "FullBand frequency segmentation" begin
        # one gain per scan shared across both channels
        pr = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1); freqseg = FullBand()))
        obs = Comrade.ObservedArrayPrior(pr, arrmf)
        cs = Comrade.sitechains(obs.sitemap)
        @test all(chs -> length(chs) == 1, values(cs))
        # the parameter count equals the single-channel case (channel 2 mirrors channel 1)
        pr1 = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        @test length(obs) == length(Comrade.ObservedArrayPrior(pr1, arr1))

        # each channel-2 datum maps to the same shared gain as its channel-1 mirror
        bm = Comrade._construct_baselinemap(arrmf, obs.sitemap)
        n1 = length(Comrade.datatable(arr1))
        @test bm.indices_1[1:n1] == bm.indices_1[(n1 + 1):end]
        @test bm.indices_2[1:n1] == bm.indices_2[(n1 + 1):end]

        # Gauss-Markov FullBand: one chain per site, same density as the single-channel prior
        prgm = ArrayPrior(
            GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0); freqseg = FullBand())
        )
        obsgm = Comrade.ObservedArrayPrior(prgm, arrmf)
        @test all(chs -> length(chs) == 1, values(Comrade.sitechains(obsgm.sitemap)))
        prgm1 = ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0)))
        obsgm1 = Comrade.ObservedArrayPrior(prgm1, arr1)
        @test length(obsgm) == length(obsgm1)
        x = rand(obsgm1)
        @test logpdf(obsgm.dists, parent(x)) ≈ logpdf(obsgm1.dists, parent(x))

        # mixed segmentation across sites: AA shares the band, everyone else is per channel
        prmix = ArrayPrior(
            IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1));
            AA = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1); freqseg = FullBand())
        )
        obsmix = Comrade.ObservedArrayPrior(prmix, arrmf)
        csm = Comrade.sitechains(obsmix.sitemap)
        @test length(csm.AA) == 1
        @test all(s -> s == :AA || length(csm[s]) == 2, keys(csm))
        bmix = Comrade._construct_baselinemap(arrmf, obsmix.sitemap)
        @test length(unique(vcat(bmix.indices_1, bmix.indices_2))) == length(obsmix)

        # SEFD referencing with mixed segmentation: the band-spanning AA parameters
        # compete inside the per-channel reference cells instead of forming their own
        # always-fixed group, so the anchor count is bounded by the number of cells
        smmix = obsmix.sitemap
        fmix, _ = Comrade.reference_indices(arrmf, smmix, SEFDReference(0.0))
        @test length(fmix) ≤ length(smmix.taxis) * length(smmix.faxis)
        # uniform segmentation is unchanged: one anchor per distinct (time, channel)
        pru = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)); refant = SEFDReference(0.0))
        smu = Comrade.ObservedArrayPrior(pru, arrmf).sitemap
        fu, _ = Comrade.reference_indices(arrmf, smu, SEFDReference(0.0))
        @test length(fu) == length(unique(collect(zip(smu.tcode, smu.fcode))))
    end

    @testset "frequency-dependent param_map" begin
        ν0 = 230.0e9
        G = SingleStokesGain() do x, p
            exp(x.lg + x.α * log(p.Fr / ν0))
        end
        intprior = (
            lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
            α = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 1.0))),
        )
        intm = InstrumentModel(G, intprior)
        ointm, printm = Comrade.set_array(intm, arrmf)
        x = rand(printm)
        vis = Comrade.measurement(dvismf)
        vout = Comrade.apply_instrument(vis, ointm, (; instrument = x))

        # hand computation: p.Fr must be the DATUM's frequency, so a channel-2 datum uses
        # the channel-2 gain params and the 345 GHz spectral factor
        # independent per-point matching against the lookup's segment edges
        lgsm = Comrade.ObservedArrayPrior(intprior.lg, arrmf).sitemap
        αsm = Comrade.ObservedArrayPrior(intprior.α, arrmf).sitemap
        findval(v, sl, t, f, s) =
            parent(v)[findfirst(i -> inseg(sl, i, t, f, s), eachindex(parent(v)))]
        T = arrmf[:Ti]
        F = arrmf[:Fr]
        bl = arrmf[:sites]
        for i in eachindex(T, F, bl)
            s1, s2 = bl[i]
            g1 = exp(findval(x.lg, lgsm, T[i], F[i], s1) + findval(x.α, αsm, T[i], F[i], s1) * log(F[i] / ν0))
            g2 = exp(findval(x.lg, lgsm, T[i], F[i], s2) + findval(x.α, αsm, T[i], F[i], s2) * log(F[i] / ν0))
            @test vout[i] ≈ g1 * vis[i] * conj(g2)
        end

        # the hot path is inferred and allocation-free
        xint = map(parent ∘ Comrade.getparams, x)
        @inferred Comrade.apply_jones(vis[1], 1, ointm, xint)
        vbuf = parent(Comrade.intout(vis))
        vp = parent(vis)
        applyall!() = Comrade._apply_instrument!(vbuf, vp, ointm, xint)
        applyall!()
        @test (@allocated applyall!()) == 0
    end

    @testset "baselinemap matches brute force (ragged)" begin
        pr = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        obs = Comrade.ObservedArrayPrior(pr, arrmfr)
        smap = obs.sitemap
        bm = Comrade._construct_baselinemap(arrmfr, smap)

        # the old O(nvis * npoints) reference implementation
        T = arrmfr[:Ti]
        F = arrmfr[:Fr]
        bl = arrmfr[:sites]
        for i in eachindex(T, F, bl)
            t = T[i]
            f = F[i]
            s1, s2 = bl[i]
            i1 = findall(j -> inseg(smap, j, t, f, s1), eachindex(smap.tcode))
            i2 = findall(j -> inseg(smap, j, t, f, s2), eachindex(smap.tcode))
            @test length(i1) == 1 && bm.indices_1[i] == only(i1)
            @test length(i2) == 1 && bm.indices_2[i] == only(i2)
        end
    end
end
