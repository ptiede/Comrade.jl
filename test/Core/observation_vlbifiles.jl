using VLBIFiles
using Pyehtim
using Comrade
using Test
using Statistics


@testset "ecef_to_geodetic" begin
    # ALMA reference: ECEF ~ (2225061, -5440057, -2481681) m
    # → geodetic latitude -23.029°, longitude -67.755°, altitude ~5070 m.
    lat, lon, h = Comrade.ecef_to_geodetic(2225061.16435, -5440057.36995, -2481681.15032)
    @test rad2deg(lat) ≈ -23.0292 atol = 1.0e-3
    @test rad2deg(lon) ≈ -67.7546 atol = 1.0e-3
    @test h ≈ 5070 atol = 50.0
end


@testset "elevation_parallactic" begin
    # SMT (Mt. Graham) ECEF approx (-1828796, -5054406, 3427865); a known
    # observation of M87 (RA=187.706°, Dec=12.391°) at MJD 57854 should be
    # above the horizon for part of the day.
    xyz = (-1.8287962e6, -5.0544068e6, 3.4278652e6)
    ra_deg = 187.7059307575226
    dec_deg = 12.39112323919932
    jd = 2400000.5 + 57854 + 4.0 / 24
    el, par = Comrade.elevation_parallactic(xyz, ra_deg, dec_deg, jd)
    @test isfinite(el)
    @test isfinite(par)
    @test -π / 2 ≤ el ≤ π / 2
    @test -π ≤ par ≤ π
end


@testset "VLBIFiles vs Pyehtim — extract_vis" begin
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    obs = Pyehtim.load_uvfits_and_array(path)

    v_jl = extract_table(uvd, Visibilities(; time_average = VLBI.GapBasedScans()))
    Pyehtim.scan_average(obs)
    v_py = Comrade.extract_vis(Pyehtim.scan_average(obs))

    cj = arrayconfig(v_jl)
    cp = arrayconfig(v_py)

    # exact metadata
    @test cj.source == cp.source
    @test cj.mjd == cp.mjd
    @test cj.bandwidth ≈ cp.bandwidth rtol = 1.0e-6
    @test sites(v_jl) == sites(v_py)
    @test cj.ra ≈ cp.ra rtol = 1.0e-9
    @test cj.dec ≈ cp.dec rtol = 1.0e-9

    # match jl <-> py rows by (sites, time-rounded) and compare geometry
    kj = Dict{Tuple{Tuple{Symbol, Symbol}, Float64}, Int}()
    for i in 1:length(v_jl)
        d = v_jl[i]
        kj[(d.baseline.sites, round(d.baseline.Ti; digits = 3))] = i
    end
    matched = 0
    elev_diff = Float64[]
    par_diff = Float64[]
    for i in 1:length(v_py)
        d = v_py[i]
        j = get(kj, (d.baseline.sites, round(d.baseline.Ti; digits = 3)), 0)
        j == 0 && continue
        matched += 1
        push!(elev_diff, abs(v_jl[j].baseline.elevation[1] - d.baseline.elevation[1]))
        push!(elev_diff, abs(v_jl[j].baseline.elevation[2] - d.baseline.elevation[2]))
        push!(par_diff, abs(v_jl[j].baseline.parallactic[1] - d.baseline.parallactic[1]))
        push!(par_diff, abs(v_jl[j].baseline.parallactic[2] - d.baseline.parallactic[2]))
    end
    @test matched > length(v_py) ÷ 2

    # Geometry tolerances: ehtim and our pure-Julia helper both use mean
    # GMST + the standard parallactic atan2 formula but differ in their
    # sidereal-time numerics, so a few-arcminute residual is expected.
    @test maximum(elev_diff) < 0.01    # rad ≈ 0.6°
    @test maximum(par_diff) < 0.02     # rad ≈ 1.1°
    @test mean(elev_diff) < 0.005
    @test mean(par_diff) < 0.005
end


@testset "VLBIFiles extract_amp" begin
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    a = extract_table(uvd, VisibilityAmplitudes(; time_average = VLBI.GapBasedScans()))
    @test length(a) > 0
    @test all(isfinite, a[:measurement])
    @test all(a[:measurement] .≥ 0)
    @test all(a[:noise] .> 0)
end


@testset "VLBIFiles extract_coherency" begin
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    c = extract_table(uvd, Coherencies(; time_average = VLBI.GapBasedScans()))
    @test length(c) > 0
    @test eltype(c[:measurement]) <: AbstractMatrix{<:Complex}
    @test all(d -> d.polbasis == (CirBasis(), CirBasis()), c[:baseline])
end


@testset "VLBIFiles extract_cphase + extract_lcamp" begin
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    cp = extract_table(uvd, ClosurePhases(; time_average = VLBI.GapBasedScans()))
    lc = extract_table(uvd, LogClosureAmplitudes(; time_average = VLBI.GapBasedScans()))
    @test length(cp) > 0
    @test length(lc) > 0
    @test all(isfinite, cp[:measurement])
    @test all(isfinite, lc[:measurement])
end


@testset "load_array_txt + reset_mounts!" begin
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)

    # ehtim-style array.txt parsing
    tmp = tempname()
    try
        open(tmp, "w") do io
            println(io, "# site X Y Z SEFDR SEFDL DRre DRim DLre DLim FR_PAR FR_ELEV FR_OFF")
            println(io, "AA  0 0 0  9999  8888  0 0  0 0  1.0  0.5  -30.0")
        end
        loaded = Comrade.load_array_txt(tmp)
        @test loaded[:AA].SEFD1 == 9999.0
        @test loaded[:AA].SEFD2 == 8888.0
        @test loaded[:AA].fr_parallactic == 1.0
        @test loaded[:AA].fr_elevation == 0.5
        @test loaded[:AA].fr_offset ≈ deg2rad(-30.0)

        # apply via reset_mounts!
        v = extract_table(uvd, Visibilities(; time_average = VLBI.GapBasedScans()))
        reset_mounts!(v, loaded)
        i = findfirst(==(:AA), arrayconfig(v).tarr.sites)
        @test arrayconfig(v).tarr.SEFD1[i] == 9999.0
        @test arrayconfig(v).tarr.fr_offset[i] ≈ deg2rad(-30.0)
    finally
        rm(tmp; force = true)
    end
end


@testset "extract_table dispatch on UVData" begin
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    v1 = extract_table(uvd, Visibilities(; time_average = VLBI.GapBasedScans()))
    v2 = extract_table(uvd, Visibilities(; time_average = VLBI.GapBasedScans()))
    @test length(v1) == length(v2)
    @test v1[:measurement] ≈ v2[:measurement]
end


@testset "time_average / frequency_average kwargs" begin
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    v_avg = extract_table(uvd, Visibilities(; time_average = VLBI.GapBasedScans()))
    @test length(v_avg) > 0

    a_avg = extract_table(uvd, VisibilityAmplitudes(; time_average = VLBI.GapBasedScans()))
    @test length(a_avg) == length(v_avg)

    cp_avg = extract_table(uvd, ClosurePhases(; time_average = VLBI.GapBasedScans()))
    @test length(cp_avg) > 0
end
