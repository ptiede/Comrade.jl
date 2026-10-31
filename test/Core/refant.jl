using Distributions

@testset "Referencing schemes" begin

    # A three-scan schedule whose first and last scans share only GL, so the sites
    # split into an early and a late group. Sites appear a different number of times,
    # which is what separates the per-stamp and per-track pin counts below.
    tsegs = [
        (0.0, [:AA, :AX, :GL]),
        (1.0, [:AA, :AX, :GL]),
        (2.0, [:GL, :LM, :SW]),
    ]
    times = Comrade.IntegrationTime{Int, Float64}[]
    sites = Symbol[]
    for (t, ss) in tsegs, s in ss
        push!(times, Comrade.IntegrationTime(1, t, 0.1))
        push!(sites, s)
    end
    freqs = fill(Comrade.FrequencyChannel(230.0e9, 1.0e9, 1), length(sites))
    scan_lookup = Comrade.SiteLookup(times, freqs, sites)

    # the same array with one entry per site, as a track segmentation produces
    trk_sites = unique(sites)
    trk_lookup = Comrade.SiteLookup(
        [Comrade.IntegrationTime(1, 1.0, 2.0) for _ in trk_sites],
        [Comrade.FrequencyChannel(230.0e9, 1.0e9, 1) for _ in trk_sites],
        trk_sites,
    )

    @testset "construction" begin
        r = MultiReference([:AA, :GL], 0.0)
        @test r.sites == [:AA, :GL]
        @test r.value == 0.0
        # any iterable of names is accepted and normalized to symbols
        @test MultiReference((:AA, :GL), 0.0).sites == [:AA, :GL]
        @test MultiReference(["AA", "GL"], 0.0).sites == [:AA, :GL]
        rp = MultiReference{Vector{Symbol}}([:AA, :GL], 0.0)
        @test rp.sites == r.sites && rp.value == r.value
        @test MultiReference([:AA], 0) isa MultiReference{Vector{Symbol}, Int}

        @test_throws "at least one site" MultiReference(Symbol[], 0.0)
        @test_throws "must be distinct" MultiReference([:AA, :AA], 0.0)
    end

    @testset "track segmentation pins one entry per site" begin
        inds, vals = Comrade.reference_indices(
            nothing, trk_lookup, MultiReference([:AA, :LM], 0.0)
        )
        @test trk_lookup.sites[inds] == [:AA, :LM]
        @test vals == [0.0, 0.0]
        @test issorted(inds)

        # one pin per site is what separates this from SingleReference
        sinds, _ = Comrade.reference_indices(
            nothing, trk_lookup, SingleReference(:AA, 0.0)
        )
        @test length(sinds) == 1
        @test length(inds) == 2
    end

    @testset "scan segmentation pins every stamp of each site" begin
        inds, vals = Comrade.reference_indices(
            nothing, scan_lookup, MultiReference([:GL, :SW], 1.5)
        )
        @test sort(scan_lookup.sites[inds]) == [:GL, :GL, :GL, :SW]
        @test all(==(1.5), vals)
    end

    @testset "build_dist folds the pins into the fixed indices" begin
        sitedists = NamedTuple(
            s => IIDSitePrior(TrackSeg(), Normal(0.0, 1.0)) for s in trk_sites
        )
        d = Comrade.build_dist(
            sitedists, trk_lookup, nothing, MultiReference([:AA, :GL], 0.0), nothing
        )
        @test d isa Comrade.PartiallyConditionedDist
        @test length(d) == length(trk_sites)
        @test sort(trk_lookup.sites[d.fixed_index]) == [:AA, :GL]
        @test length(d.variate_index) == length(trk_sites) - 2
        x = rand(d)
        @test x[d.fixed_index] == [0.0, 0.0]
        @test isfinite(logpdf(d, x))
    end

    @testset "a first-stamp pin must agree with the reference value" begin
        # AA is both reference-fixed and init-pinned here; disagreeing values would
        # silently change which gauge the model is in.
        agree = NamedTuple(
            s => IIDSitePrior(
                    TrackSeg(), Normal(0.0, 1.0); init = s === :AA ? FixedInit(0.0) : nothing
                ) for s in trk_sites
        )
        d = Comrade.build_dist(
            agree, trk_lookup, nothing, MultiReference([:AA, :GL], 0.0), nothing
        )
        @test sort(trk_lookup.sites[d.fixed_index]) == [:AA, :GL]

        clash = NamedTuple(
            s => IIDSitePrior(
                    TrackSeg(), Normal(0.0, 1.0); init = s === :AA ? FixedInit(0.7) : nothing
                ) for s in trk_sites
        )
        @test_throws "conflicts with the reference value" Comrade.build_dist(
            clash, trk_lookup, nothing, MultiReference([:AA, :GL], 0.0), nothing
        )
    end

    @testset "an unknown site is an error, not a silent no-op" begin
        @test_throws "not among this parameter's sites" Comrade.reference_indices(
            nothing, trk_lookup, MultiReference([:AA, :XX], 0.0)
        )
    end

end
