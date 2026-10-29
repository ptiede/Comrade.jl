using VLBIFiles
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

        t = VisTable(uvd, IPol())
        reset_mounts!(t, loaded)
        st = info(t).stations
        i = findfirst(==(:AA), st.name)
        @test st.mount[i].parallactic == 1.0 && st.mount[i].elevation == 0.5
        @test all(≈(deg2rad(-30.0)), st.feed_angles[i])
    finally
        rm(tmp; force = true)
    end
end
