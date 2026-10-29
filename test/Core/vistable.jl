using Comrade, Test, LinearAlgebra, Random, Distributions
using VLBIFiles

function _vistable_info(; frequency = [227.0e9, 229.0e9])
    return ArrayInfo(;
        stations = (;
            name = [:AA, :AP, :SM, :JC, :LM],
            position = [
                (2225061.0, -5440057.0, -2481681.0), (2225039.0, -5441197.0, -2479303.0),
                (-5464584.0, -2493001.0, 2150654.0), (-5464555.0, -2492928.0, 2150797.0),
                (-768713.0, -5988541.0, 2063275.0),
            ],
            mount = [MountNasmythL(), MountAltAz(), MountAltAz(), MountAltAz(), MountNasmythR()],
            feeds = fill((RPol(), LPol()), 5),
            feed_angles = fill((0.0, 0.0), 5),
        ),
        channels = (; frequency, width = fill(2.0e9, length(frequency))),
        scans = (; start = [0.9], stop = [3.1]),
        source = Source("M87", 187.70593, 12.39112),
        epoch = Comrade.Date(2017, 4, 11),
    )
end

@testset "VisTable" begin
    ainfo = _vistable_info()
    n = 6
    Ti = [2.0, 1.0, 1.0, 2.0, 1.0, 2.0]
    ant1 = [:AA, :AA, :AP, :AP, :AA, :AA]
    ant2 = [:AP, :SM, :SM, :SM, :LM, :SM]
    rng = Random.Xoshiro(1)
    V = complex.(randn(rng, n, 2), randn(rng, n, 2))
    σ = fill(0.1, n, 2)
    σ[3, 2] = Inf
    V[4, 1] = NaN
    u, v = 1.0e3 .* randn(rng, n), 1.0e3 .* randn(rng, n)
    t = VisTable(IPol(), V, σ; u, v, Ti, ant1, ant2, info = ainfo)

    @testset "construction" begin
        d = domain(t)
        @test d.Ti == sort(Ti)
        @test issorted(collect(zip(d.Ti, d.ant1, d.ant2)))
        @test count(iszero, parent(weight(t))) == 2
        @test all(isfinite, parent(measurement(t)))
        @test isinf(maximum(parent(noise(t))))
        @test size(d.U) == (n, 2)
        @test d.U[:, 2] ≈ ComradeBase.coords(d).u .* 229.0e9 ./ ComradeBase.speed_of_light
        @test baselines(t)[1] == (:AA, :SM)
        @test sites(t) == [:AA, :AP, :SM, :LM]
        @test size(stations(t).elevation) == (5, 2)
        @test all(x -> -π / 2 <= x <= π / 2, stations(t).elevation)
        @test info(t) === ainfo
    end

    @testset "construction errors" begin
        @test_throws "give either the noise" VisTable(IPol(), V; u, v, Ti, ant1, ant2, info = ainfo)
        @test_throws "joins a station to itself" VisTable(IPol(), V, σ; u, v, Ti, ant1 = fill(:AA, n), ant2 = fill(:AA, n), info = ainfo)
        @test_throws "is not in the array" VisTable(IPol(), V, σ; u, v, Ti, ant1 = fill(:XX, n), ant2, info = ainfo)
        @test_throws DimensionMismatch VisTable(IPol(), V[:, 1], σ[:, 1]; u, v, Ti, ant1, ant2, info = ainfo)
        @test_throws "weights must be finite" VisTable(IPol(), V; weight = -σ, u, v, Ti, ant1, ant2, info = ainfo)
        @test_throws "missing" ArrayInfo(; stations = (; name = [:AA]), channels = ainfo.channels, scans = ainfo.scans, source = ainfo.source, epoch = ainfo.epoch)
        @test_throws "strictly increasing" _vistable_info(; frequency = [229.0e9, 227.0e9])
        st = ainfo.stations
        stationtable(; kw...) = merge((; name = st.name, position = st.position, mount = st.mount, feeds = st.feeds, feed_angles = st.feed_angles), (; kw...))
        @test ainfo.tracks.start == [0.9] && ainfo.tracks.stop == [3.1]
        tracks_with(tracks) = ArrayInfo(; ainfo.stations, ainfo.channels, scans = (; start = [0.9, 5.0], stop = [3.1, 6.0]), ainfo.source, ainfo.epoch, tracks)
        @test length(tracks_with((; start = [0.0, 4.0], stop = [3.5, 7.0])).tracks) == 2
        @test_throws "in order and disjoint" tracks_with((; start = [0.0, 3.0], stop = [4.0, 7.0]))
        @test_throws "is not inside a track" tracks_with((; start = [0.0], stop = [4.0]))
        info_with(stations) = ArrayInfo(; stations, channels = ainfo.channels, scans = ainfo.scans, source = ainfo.source, epoch = ainfo.epoch)
        @test_throws "feeds and" info_with(stationtable(; feed_angles = fill((0.0,), 5)))
        @test_throws "feeds are polarizations" info_with(stationtable(; feeds = fill((:R, :L), 5)))
        @test_throws "AbstractMount" info_with(stationtable(; mount = fill(1.0, 5)))
        @test info_with(stationtable(; diameter = fill(12.0, 5))).stations.diameter[1] == 12.0
    end

    @testset "polarization and amplitudes" begin
        @test_throws MethodError VisTable(:I, V, σ; u, v, Ti, ant1, ant2, info = ainfo)
        @test_throws "the feeds 1 and 2" Correlation(1, 3)
        t12 = VisTable(Correlation(1, 2), V, σ; u, v, Ti, ant1, ant2, info = ainfo)
        @test t12.pol === Correlation(1, 2)
        a = amplitudes(t)
        @test parent(measurement(a)) == abs.(parent(measurement(t)))
        @test weight(a) == weight(t)
        ad = amplitudes(t; debias = true)
        w = parent(weight(t))
        @test parent(measurement(ad))[w .> 0] ≈ sqrt.(max.(abs2.(parent(measurement(t))[w .> 0]) .- inv.(w[w .> 0]), 0))
        @test_throws "then take their amplitudes" average(a, ByScan())
        @test_throws "scalar complex visibilities" closure_phases(a)
        @test_throws "scalar complex visibilities" amplitudes(a)
    end

    @testset "coherencies" begin
        info1 = _vistable_info(; frequency = [230.0e9])
        C = complex.(randn(rng, n, 1, 2, 2), randn(rng, n, 1, 2, 2))
        σc = fill(0.1, n, 1, 2, 2)
        σc[1, 1, 1, 2] = NaN
        tc = VisTable(Correlations(), C, σc; u, v, Ti, ant1, ant2, info = info1)
        @test eltype(measurement(tc)) <: Comrade.StaticArraysCore.SMatrix{2, 2, ComplexF64}
        @test count(iszero, parent(parent(weight(tc)))) == 1
        @test size(drop(tc)) == size(tc)
        tca = add_fractional_noise(tc, 0.1)
        @test all(parent(parent(weight(tca))) .<= parent(parent(weight(tc))))
        @test_throws "scalar complex visibilities" closure_phases(tc)
        @test_throws "scalar complex visibilities" amplitudes(tc)
    end

    @testset "operations" begin
        t1 = t[Fr = At(229.0e9)]
        @test size(t1) == (n, 1)
        @test domain(t1).U ≈ domain(t).U[:, 2:2]
        t2 = t[Node = 2:4]
        @test size(t2) == (3, 2)
        @test domain(t2).Ti == domain(t).Ti[2:4]

        tf = flag(r -> r.baseline == (:AA, :SM), t)
        flagged = baselines(t) .== Ref((:AA, :SM))
        @test all(iszero, parent(weight(tf))[flagged, :])
        @test parent(weight(tf))[.!flagged, :] == parent(weight(t))[.!flagged, :]
        tfd = drop(flag(r -> r.Ti < 1.5, t))
        @test all(>(1.5), domain(tfd).Ti)
        @test size(tfd, 1) == count(>(1.5), Ti)

        ta = add_fractional_noise(t, 0.1)
        w0, wa = parent(weight(t)), parent(weight(ta))
        @test iszero.(wa) == iszero.(w0)
        @test inv(wa[1, 1]) ≈ inv(w0[1, 1]) + 0.01 * abs2(parent(measurement(t))[1, 1])
        @test beamsize(t) ≈ inv(maximum(hypot.(domain(t).U, domain(t).V)[parent(weight(t)) .> 0]))
    end

    @testset "joining" begin
        a = info(t)
        info2 = ArrayInfo(;
            stations = (;
                name = [:AP, :AA, :PV],
                position = [a.stations.position[2], a.stations.position[1], (5088968.0, -301681.0, 3825012.0)],
                mount = [MountAltAz(), MountNasmythL(), MountNasmythL()],
                feeds = fill((RPol(), LPol()), 3), feed_angles = fill((0.0, 0.0), 3),
            ),
            channels = (; frequency = [229.0e9, 231.0e9], width = fill(2.0e9, 2)),
            scans = (; start = [0.5, 1.0], stop = [1.15, 1.5]),
            source = a.source, epoch = Comrade.Date(2017, 4, 12),
        )
        V2 = complex.(randn(rng, 3, 2), randn(rng, 3, 2))
        t2 = VisTable(
            IPol(), V2, fill(0.2, 3, 2); u = 1.0e3 .* randn(rng, 3), v = 1.0e3 .* randn(rng, 3),
            Ti = [1.0, 1.1, 1.2], ant1 = [:AA, :AP, :AA], ant2 = [:PV, :PV, :AP], info = info2,
        )
        j = vcat(t, t2)
        ja = info(j)
        @test ja.epoch == a.epoch
        @test ja.stations.name == [:AA, :AP, :SM, :JC, :LM, :PV]
        @test ja.channels.frequency == [227.0e9, 229.0e9, 231.0e9]
        @test ja.scans.start ≈ [0.9, 24.5] && ja.scans.stop ≈ [3.1, 25.5]
        @test ja.tracks.start ≈ [0.9, 24.5] && ja.tracks.stop ≈ [3.1, 25.5]
        @test size(j) == (n + 3, 3)
        @test domain(j).Ti ≈ vcat(domain(t).Ti, [25.0, 25.1, 25.2])
        w, m = parent(weight(j)), parent(measurement(j))
        old, new = 1:n, (n + 1):(n + 3)
        @test w[old, 1:2] == parent(weight(t)) && m[old, 1:2] == parent(measurement(t))
        @test w[new, 2:3] == parent(weight(t2)) && m[new, 2:3] == parent(measurement(t2))
        @test all(iszero, w[old, 3]) && all(iszero, w[new, 1])
        @test baselines(j) == vcat(baselines(t), baselines(t2))
        @test domain(j).U[new, 2:3] ≈ domain(t2).U
        el = stations(j).elevation
        @test el[Ant = At(:PV), Ts = At(domain(j).Ti[n + 1])] ≈ stations(t2).elevation[Ant = At(:PV), Ts = At(1.0)]
        @test parent(el)[1:5, 1:2] ≈ parent(stations(t).elevation)

        @test_throws "measure the same channel" vcat(t, t)
        split = vcat(t[Fr = At(227.0e9)], t[Fr = At(229.0e9)])
        @test size(split) == size(t)
        @test parent(weight(split)) == parent(weight(t)) && parent(measurement(split)) == parent(measurement(t))
        @test domain(split).U ≈ domain(t).U
        lo1 = VisTable(IPol(), V[:, 1], σ[:, 1]; u, v, Ti, ant1, ant2, info = _vistable_info(; frequency = [227.0e9]))
        hi1(δ) = VisTable(IPol(), V[:, 2], σ[:, 2]; u = u .* (1 + δ), v, Ti, ant1, ant2, info = _vistable_info(; frequency = [229.0e9]))
        @test domain(vcat(lo1, hi1(4.0e-4))).U[:, 1] ≈ domain(lo1).U[:, 1] .* (1 + 2.0e-4) rtol = 1.0e-6
        @test length(info(vcat(lo1, hi1(4.0e-4))).tracks) == 1
        @test_throws "u, v differing by more than" vcat(lo1, hi1(1.0e-2))
        @test_throws "amplitude tables cannot be joined" vcat(amplitudes(t), amplitudes(t2))
        @test_throws "cannot be joined" vcat(t, VisTable(Correlation(1, 1), V, σ; u, v, Ti, ant1, ant2, info = a))
        moved = ArrayInfo(; stations = merge(Comrade.Tables.columntable(a.stations), (; position = reverse(a.stations.position))), a.channels, a.scans, a.source, a.epoch)
        @test_throws "station AA differs" vcat(t, VisTable(IPol(), V, σ; u, v, Ti, ant1, ant2, info = moved))
        near = ArrayInfo(; a.stations, channels = (; frequency = [227.5e9], width = [2.0e9]), a.scans, a.source, a.epoch)
        @test_throws "overlap" vcat(t, VisTable(IPol(), V[:, 1], σ[:, 1]; u, v, Ti, ant1, ant2, info = near))
    end

    @testset "averaging" begin
        d = domain(t)
        ta = average(t, _ -> 1)
        @test size(ta, 1) == length(unique(baselines(t)))
        @test times(ta) == [(minimum(Ti) + maximum(Ti)) / 2]
        w, x = parent(weight(t)), parent(measurement(t))
        sel = baselines(t) .== Ref((:AA, :SM))
        i = findfirst(==((:AA, :SM)), baselines(ta))
        @test parent(weight(ta))[i, :] ≈ vec(sum(w[sel, :]; dims = 1))
        @test parent(measurement(ta))[i, :] ≈ vec(sum(w[sel, :] .* x[sel, :]; dims = 1) ./ sum(w[sel, :]; dims = 1))
        @test ComradeBase.coords(domain(ta)).u[i] ≈ sum(ComradeBase.coords(d).u[sel]) / count(sel)

        tf = average(t, ChannelBins(2))
        @test size(tf) == (n, 1)
        @test info(tf).channels.frequency == [228.0e9]
        @test info(tf).channels.width == [4.0e9]
        @test parent(weight(tf))[:, 1] ≈ vec(sum(w; dims = 2))
        @test parent(measurement(average(average(t, ChannelBins(2)), ByScan()))) ≈ parent(measurement(average(t, ByScan(), ChannelBins(2))))
        @test parent(measurement(average(t, ByChannel(t -> [2, 1])))) ≈ parent(measurement(t))
        @test_throws "run length must divide" average(t, ChannelBins(3))
        @test_throws "positive number of channels" ChannelBins(0)
        late = VisTable(IPol(), V, σ; u, v, Ti = Ti .+ 5, ant1, ant2, info = ainfo)
        @test_throws "outside every scan" average(late, ByScan())
        @test_throws DimensionMismatch average(t, _ -> [1, 2, 3])
    end

    @testset "closures" begin
        names = ainfo.stations.name
        full = [(a, b) for (i, a) in pairs(names) for (j, b) in pairs(names) if j > i]
        square = [(:AA, :AP), (:AP, :SM), (:SM, :JC), (:JC, :AA)]
        ring = [(:AA, :AP), (:AP, :SM), (:SM, :JC), (:JC, :LM), (:LM, :AA)]
        spur = [(:AP, :AA), (:SM, :AP), (:AA, :SM), (:JC, :SM)]
        groups = [(1.0, full), (2.0, square), (2.5, ring), (3.0, spur)]
        cTi = [t for (t, prs) in groups for _ in prs]
        c1 = [a for (_, prs) in groups for (a, _) in prs]
        c2 = [b for (_, prs) in groups for (_, b) in prs]
        m = length(cTi)
        info1 = _vistable_info(; frequency = [230.0e9])
        Vc = complex.(1 .+ rand(rng, m), randn(rng, m))
        tt = VisTable(IPol(), Vc, fill(0.01, m); u = randn(rng, m), v = randn(rng, m), Ti = cTi, ant1 = c1, ant2 = c2, info = info1)
        cp, la = closure_phases(tt), log_closure_amplitudes(tt)
        perstamp(ct) = [count(==(t), domain(ct).Ti) for (t, _) in groups]
        # E − V + C cycles; E − V gain-free amplitude ratios, one more without an odd cycle
        @test perstamp(cp) == [10 - 5 + 1, 4 - 4 + 1, 5 - 5 + 1, 4 - 4 + 1]
        @test perstamp(la) == [10 - 5, 4 - 4 + 1, 5 - 5, 4 - 4]
        @test rank(Matrix(getfield(cp, :design))) == length(cp)
        @test rank(Matrix(getfield(la, :design))) == length(la)
        @test all(s -> length(s) == 4, domain(cp).stations[domain(cp).Ti .== 2.0])

        g = Dict(s => exp(complex(randn(rng), randn(rng))) for s in names)
        d = domain(tt)
        x = vec(parent(measurement(tt)))
        xg = [x[i] * g[names[d.ant1[i]]] * conj(g[names[d.ant2[i]]]) for i in eachindex(x)]
        wrap(y) = mod2pi(y + π) - π
        @test all(abs.(wrap.(getfield(cp, :design) * angle.(xg) .- parent(measurement(cp)))) .< 1.0e-10)
        @test getfield(la, :design) * log.(abs.(xg)) ≈ parent(measurement(la))
        @test noise(cp) ≈ sqrt.(diag(getfield(cp, :covariance)))
    end

    @testset "VLBIFiles" begin
        file = joinpath(@__DIR__, "..", "..", "examples", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits")
        uvd = VLBIFiles.load(VLBIFiles.UVData, file)
        raw = VisTable(uvd, IPol())
        @test_throws MethodError VisTable(uvd, :I)
        rr = VisTable(uvd, Correlation(1, 1))
        ll = VisTable(uvd, Correlation(2, 2))
        @test size(rr) == size(ll) == size(raw)
        wr, wl = parent(weight(rr)), parent(weight(ll))
        both = (wr .> 0) .& (wl .> 0)
        @test any(both)
        @test parent(measurement(raw))[both] ≈ ((parent(measurement(rr)) .+ parent(measurement(ll))) ./ 2)[both]
        @test parent(weight(raw))[both] ≈ (4 ./ (inv.(wr) .+ inv.(wl)))[both]
        hi = VisTable(VLBIFiles.load(VLBIFiles.UVData, replace(file, "_lo_" => "_hi_")), IPol())
        lohi = vcat(raw, hi)
        records(t) = collect(zip(domain(t).Ti, baselines(t)))
        @test size(lohi) == (length(union(records(raw), records(hi))), 2)
        @test count(>(0), parent(weight(lohi))) == count(>(0), parent(weight(raw))) + count(>(0), parent(weight(hi)))
        @test sort(vec(parent(measurement(lohi))[parent(weight(lohi)) .> 0]); by = reim) ==
            sort(vcat(vec(parent(measurement(raw))[parent(weight(raw)) .> 0]), vec(parent(measurement(hi))[parent(weight(hi)) .> 0])); by = reim)

        t = average(raw, ByScan())
        @test 0 < count(>(0), parent(weight(t))) < count(>(0), parent(weight(raw)))
        @test all(τ -> any(k -> info(t).scans.start[k] <= τ <= info(t).scans.stop[k], eachindex(info(t).scans.start)), domain(t).Ti)
    end
end
