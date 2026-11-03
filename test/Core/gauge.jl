using Distributions
using LinearAlgebra
using Random

@testset "Phase gauge solver" begin

    GAUGE_MJD = 57849
    GAUGE_FREQ = 230.0e9

    _chan() = Comrade.FrequencyChannel(GAUGE_FREQ, 1.0e9, 1)

    # One entry per (stamp, site), as an integration-time segmentation produces. Returns the
    # lookup and the (stamp, site) => entry index map the pin lists below are written with.
    function integ_lookup(times, sites_at; dt = 0.05)
        ts = Comrade.IntegrationTime{Int, Float64}[]
        ss = Symbol[]
        idx = Dict{Tuple{Int, Symbol}, Int}()
        for τ in eachindex(times, sites_at), s in sites_at[τ]
            push!(ts, Comrade.IntegrationTime(GAUGE_MJD, times[τ], dt))
            push!(ss, s)
            idx[(τ, s)] = length(ss)
        end
        return Comrade.SiteLookup(ts, fill(_chan(), length(ss)), ss), idx
    end

    # One entry per site, covering every stamp, as a track segmentation produces.
    function track_lookup(times, sites_at)
        ss = unique(reduce(vcat, sites_at))
        lo, hi = extrema(times)
        span = 2 * (hi - lo + 1)
        ts = [Comrade.IntegrationTime(GAUGE_MJD, (lo + hi) / 2, span) for _ in ss]
        idx = Dict(s => i for (i, s) in pairs(ss))
        return Comrade.SiteLookup(ts, fill(_chan(), length(ss)), ss), idx
    end

    stamplist(times) = [(t, GAUGE_FREQ) for t in times]

    # The entry a per-stamp SEFDReference pins: the lowest-SEFD site present at each stamp.
    sefd_pins(sites_at, sefd, idx) =
        [idx[(τ, argmin(s -> sefd[s], sites_at[τ]))] for τ in eachindex(sites_at)]

    # The entry a FixedInit pins: each site's earliest stamp.
    function init_pins(sites_at, idx)
        seen = Dict{Symbol, Int}()
        for τ in eachindex(sites_at), s in sites_at[τ]
            haskey(seen, s) || (seen[s] = idx[(τ, s)])
        end
        return sort!(collect(values(seen)))
    end

    @testset "one per-stamp term is the classical reference antenna" begin
        times = [0.0, 1.0, 2.0]
        sites_at = [[:AA, :AX, :GL], [:AX, :GL, :LM], [:GL, :LM, :SW]]
        sefd = (AA = 100.0, AX = 300.0, GL = 15000.0, LM = 11000.0, SW = 12000.0)
        smap, idx = integ_lookup(times, sites_at; dt = 0.5)

        pins = Comrade.gauge_pins(
            (gp = (smap = smap, fixed = Int[]),), stamplist(times), sites_at, sefd
        )
        # one gauge freedom per stamp, and the entry taken is the one SEFDReference takes
        @test length(pins) == length(times)
        @test all(p -> p[1] === :gp, pins)
        @test sort(last.(pins)) == sort(sefd_pins(sites_at, sefd, idx))

        # once those entries are pinned nothing is left over
        @test isempty(
            Comrade.gauge_pins(
                (gp = (smap = smap, fixed = sefd_pins(sites_at, sefd, idx)),),
                stamplist(times), sites_at, sefd
            )
        )
    end

    @testset "an offset plus a residual on the M87 schedule" begin
        # The 2022 band-3 M87 track: the first twelve stamps are the three sites that see
        # the source before ALMA rises.
        schedule = [
            (21.5083, "GL,NN,PV"), (21.6417, "GL,NN,PV"), (21.9917, "GL,NN,PV"),
            (22.125, "GL,NN,PV"), (22.5917, "GL,NN,PV"), (22.725, "GL,NN,PV"),
            (23.075, "GL,NN,PV"), (23.2083, "GL,NN,PV"), (23.675, "GL,NN,PV"),
            (23.8083, "GL,NN,PV"), (24.1583, "GL,NN,PV"), (24.2917, "GL,NN,PV"),
            (24.875, "AA,AX,GL,NN,PV"), (25.0083, "AA,AX,GL,NN,PV"),
            (25.325, "AA,AX,GL,NN,PV"), (25.6417, "AA,AX,GL,LM,NN,PV"),
            (25.975, "AA,AX,GL,KT,LM,NN,PV"), (26.6417, "AA,AX,GL,KT,LM,NN,PV"),
            (27.0417, "AA,AX,GL,KT,NN,PV"), (27.175, "AA,AX,GL,KT,NN,PV"),
            (27.6764, "AA,AX,GL,KT,LM,MG,NN,PV"), (27.975, "AA,AX,GL,KT,LM,MG,NN,PV"),
            (28.325, "AA,AX,GL,KT,LM,MG,NN,PV"), (28.4583, "AA,AX,GL,KT,LM,MG,NN,PV"),
            (28.775, "AA,GL,KT,LM,MG,NN,PV"), (28.8917, "AA,GL,KT,LM,MG,NN,PV"),
            (29.225, "AA,AX,GL,KT,LM,MG,MM,NN,PV"), (29.3583, "AA,AX,GL,KT,LM,MG,MM,NN,PV"),
            (29.675, "AA,AX,GL,KT,LM,MG,MM,NN,PV,SW"),
            (29.8083, "AA,AX,GL,KT,LM,MG,MM,NN,PV,SW"),
            (30.125, "AA,GL,KT,LM,MG,MM,PV,SW"), (30.2583, "AA,AX,GL,KT,LM,MG,MM,PV"),
            (31.0083, "AA,AX,GL,KT,MG,MM,SW"), (31.1417, "AA,AX,GL,KT,MG,MM"),
            (32.2417, "AA,AX,GL,KT,LM,MG,MM,SW"), (32.375, "AA,AX,GL,KT,LM,MG,MM"),
            (34.5583, "KT,LM,MG,MM"), (35.225, "KT,LM,MG,MM,SW"),
            (35.6917, "KT,LM,MG,MM,SW"), (36.1583, "KT,LM,MG,MM,SW"),
        ]
        times = first.(schedule)
        sites_at = [Symbol.(split(s, ',')) for s in last.(schedule)]
        # GL, NN and PV are given the same SEFD so that the choice among them falls through
        # to the last tiebreak, the number of stamps a site is present at.
        sefd = (
            AA = 100.0, AX = 300.0, KT = 8000.0, MG = 9000.0, MM = 10000.0, LM = 11000.0,
            SW = 12000.0, GL = 15000.0, NN = 15000.0, PV = 15000.0,
        )

        res, ridx = integ_lookup(times, sites_at)
        off, oidx = track_lookup(times, sites_at)
        resfixed = sort!(unique!(vcat(sefd_pins(sites_at, sefd, ridx), init_pins(sites_at, ridx))))
        terms(offfixed) = (
            gp1μ = (smap = off, fixed = offfixed),
            gp1 = (smap = res, fixed = resfixed),
        )

        # ALMA is the only offset the reference scheme pins, and that leaves the common
        # level of the three sites that observe before it rises undetermined.
        pins = Comrade.gauge_pins(terms([oidx[:AA]]), stamplist(times), sites_at, sefd)
        @test length(pins) == 1
        name, i = only(pins)
        @test name === :gp1μ
        @test off.sites[i] ∈ (:GL, :NN, :PV)
        # GL is present at more stamps than NN or PV
        @test off.sites[i] === :GL

        # applying the pin identifies the parameterization
        @test isempty(
            Comrade.gauge_pins(
                terms(sort([oidx[:AA], i])), stamplist(times), sites_at, sefd
            )
        )
    end

    @testset "the coarsest term in a group carries the pin" begin
        times = [0.0]
        sites_at = [[:AA, :LM]]
        sefd = (AA = 100.0, LM = 11000.0)
        res, ridx = integ_lookup(times, sites_at)
        off, oidx = track_lookup(times, sites_at)

        # AA's residual and LM's offset are pinned, so the one remaining group holds a
        # track-long entry and an integration-long one.
        terms = (
            gp1μ = (smap = off, fixed = [oidx[:LM]]),
            gp1 = (smap = res, fixed = [ridx[(1, :AA)]]),
        )
        @test Comrade.gauge_pins(terms, stamplist(times), sites_at, sefd) ==
            [(:gp1μ, oidx[:AA])]
        # the preference is the caller's to change
        @test Comrade.gauge_pins(
            terms, stamplist(times), sites_at, sefd; prefer = e -> (e.dt,)
        ) == [(:gp1, ridx[(1, :LM)])]
    end

    @testset "a datum tying three terms at once is rejected" begin
        times = [0.0, 1.0]
        sites_at = [[:AA, :LM], [:AA, :LM]]
        sefd = (AA = 100.0, LM = 11000.0)
        off1, _ = track_lookup(times, sites_at)
        off2, _ = track_lookup(times, sites_at)
        terms = (
            gpa = (smap = off1, fixed = Int[]),
            gpb = (smap = off2, fixed = Int[]),
        )
        @test_throws "leaves 2 gauge terms free at once" Comrade.gauge_pins(
            terms, stamplist(times), sites_at, sefd
        )
    end

    @testset "input checks" begin
        times = [0.0]
        sites_at = [[:AA, :LM]]
        sefd = (AA = 100.0, LM = 11000.0)
        off, _ = track_lookup(times, sites_at)
        @test_throws "fixed indices outside" Comrade.gauge_pins(
            (gp = (smap = off, fixed = [7]),), stamplist(times), sites_at, sefd
        )
        @test_throws DimensionMismatch Comrade.gauge_pins(
            (gp = (smap = off, fixed = Int[]),), stamplist([0.0, 1.0]), sites_at, sefd
        )
    end

    @testset "explicit entry pins compose with a named scheme" begin
        times = [0.0, 1.0]
        sites_at = [[:AA, :LM], [:AA, :LM]]
        smap, idx = integ_lookup(times, sites_at; dt = 0.5)
        r = Comrade.CompositeReference(
            SingleReference(:AA, 0.0), Comrade.EntryReference([idx[(2, :LM)]], 0.0)
        )
        inds, vals = Comrade.reference_indices(nothing, smap, r)
        @test issorted(inds)
        @test sort(smap.sites[inds]) == [:AA, :AA, :LM]
        @test all(==(0.0), vals)

        # a scheme that pins nothing contributes nothing
        inds2, _ = Comrade.reference_indices(
            nothing, smap, Comrade.CompositeReference(NoReference(), SingleReference(:AA, 1.0))
        )
        @test smap.sites[inds2] == [:AA, :AA]
        @test isempty(
            first(
                Comrade.reference_indices(
                    nothing, smap, Comrade.CompositeReference(NoReference())
                )
            )
        )

        @test_throws "to both" Comrade.reference_indices(
            nothing, smap,
            Comrade.CompositeReference(
                SingleReference(:AA, 0.0), Comrade.EntryReference([idx[(1, :AA)]], 1.0)
            )
        )
        @test_throws "outside this parameter's" Comrade.reference_indices(
            nothing, smap, Comrade.EntryReference([99], 0.0)
        )
        @test_throws "must be distinct" Comrade.EntryReference([1, 1], 0.0)
        @test Comrade.reference_value(NoReference()) === nothing
        @test Comrade.reference_value(r) == 0.0
    end

    @testset "disconnected baselines at one time give one stamp per group" begin
        @test Comrade._connected_sites([(:AA, :LM), (:GL, :PV), (:LM, :SW)]) ==
            [[:AA, :LM, :SW], [:GL, :PV]]
        @test Comrade._connected_sites([(:AA, :LM), (:LM, :GL), (:GL, :AA)]) ==
            [[:AA, :LM, :GL]]

        # AA-LM and GL-PV observe together with no baseline between them, so the time
        # carries two phase constants and a per-stamp term needs a pin in each group
        times = [0.0, 0.0]
        sites_at = [[:AA, :LM], [:GL, :PV]]
        sefd = (AA = 100.0, LM = 11000.0, GL = 15000.0, PV = 1500.0)
        smap, idx = integ_lookup(times, sites_at; dt = 0.5)
        pins = Comrade.gauge_pins(
            (gp = (smap = smap, fixed = Int[]),), stamplist(times), sites_at, sefd
        )
        @test sort(last.(pins)) == sort([idx[(1, :AA)], idx[(2, :PV)]])
    end

    @testset "pin count equals the nullity of the baseline phase design" begin
        # Dense cross-check: one row per baseline, +1 on the free entries covering the
        # first site and -1 on those covering the second. Baselines join every pair of
        # sites within a stamp and none across stamps.
        function nullity(lookups, fixed, sites_at)
            ncols = [length(l[1].sites) for l in lookups]
            offs = cumsum([0; ncols[1:(end - 1)]])
            rows = Vector{Float64}[]
            for τ in eachindex(sites_at)
                ss = sites_at[τ]
                for a in eachindex(ss), b in (a + 1):lastindex(ss)
                    row = zeros(sum(ncols))
                    for (k, (_, idx)) in pairs(lookups), (s, sgn) in ((ss[a], 1.0), (ss[b], -1.0))
                        i = get(idx, (τ, s), get(idx, s, 0))
                        (i == 0 || i in fixed[k]) && continue
                        row[offs[k] + i] += sgn
                    end
                    push!(rows, row)
                end
            end
            A = reduce(vcat, permutedims.(rows))
            nfree = sum(ncols) - sum(length, fixed)
            return nfree - rank(A)
        end

        rng = Random.Xoshiro(20260919)
        allsites = [:AA, :AX, :GL, :LM, :SW, :PV]
        sefd = (AA = 100.0, AX = 300.0, GL = 15000.0, LM = 11000.0, SW = 12000.0, PV = 1500.0)
        for _ in 1:200
            # random schedule; a time with at least four sites may split into two groups
            times = Float64[]
            sites_at = Vector{Symbol}[]
            for t in 1:rand(rng, 2:6)
                ss = Random.shuffle(rng, allsites)[1:rand(rng, 2:6)]
                if length(ss) >= 4 && rand(rng) < 0.4
                    c = rand(rng, 2:(length(ss) - 2))
                    append!(times, (t, t))
                    push!(sites_at, ss[1:c], ss[(c + 1):end])
                else
                    push!(times, t)
                    push!(sites_at, ss)
                end
            end
            res, ridx = integ_lookup(times, sites_at; dt = 0.05)
            off, oidx = track_lookup(times, sites_at)
            rfixed = rand(rng, Bool) ? init_pins(sites_at, ridx) : Int[]
            rand(rng, Bool) && (rfixed = union(rfixed, sefd_pins(sites_at, sefd, ridx)))
            ofixed = rand(rng, Bool) ? [first(values(oidx))] : Int[]
            lookups = ((off, oidx), (res, ridx))
            terms(of, rf) = (gpμ = (smap = off, fixed = of), gp = (smap = res, fixed = rf))

            pins = Comrade.gauge_pins(terms(ofixed, rfixed), stamplist(times), sites_at, sefd)
            @test length(pins) == nullity(lookups, (ofixed, rfixed), sites_at)

            ofull = union(ofixed, [i for (n, i) in pins if n === :gpμ])
            rfull = union(rfixed, [i for (n, i) in pins if n === :gp])
            @test nullity(lookups, (ofull, rfull), sites_at) == 0
            @test isempty(
                Comrade.gauge_pins(terms(ofull, rfull), stamplist(times), sites_at, sefd)
            )
        end
    end

    @testset "gauge is declared on the ArrayPrior" begin
        p = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0)); gauge = :phase)
        @test p.gauge === :phase
        @test ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0))).gauge === :none
        @test_throws "must be :none or :phase" ArrayPrior(
            IIDSitePrior(ScanSeg(), Normal(0.0, 1.0)); gauge = :amplitude
        )
    end

end
