using Dates

const ARRAYFILE = joinpath(@__DIR__, "..", "examples", "Data", "array.txt")

# The scan-averaged Stokes I visibilities of `test_data.uvfits`, with the stations' mounts
# from the example `array.txt`.
function load_test_table()
    uvd = VLBIFiles.load(VLBIFiles.UVData, joinpath(@__DIR__, "test_data.uvfits"))
    return average(reset_mounts(VisTable(uvd, IPol()), ARRAYFILE), ByScan())
end

# Noise-free scan-averaged data of two Gaussians (FWHMs 40 × 20 and 20 × 10 μas) on the
# `test_data.uvfits` coverage, and the coherencies of the polarized example file.
function load_data()
    t = load_test_table()
    fwhm = 2 * sqrt(2 * log(2))
    θ = (
        f1 = 1.0, σ1 = μas2rad(40.0) / fwhm, τ1 = 0.5, ξ1 = π / 3,
        f2 = 0.5, σ2 = μas2rad(20.0) / fwhm, τ2 = 0.5, ξ2 = π / 6,
        x = μas2rad(30.0), y = μas2rad(30.0),
    )
    skym = SkyModel(test_model, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 64, 64))
    vt = only(simulate_observation(VLBIPosterior(t; sky = skym), (; sky = θ); add_thermal_noise = false))

    amp = amplitudes(vt)
    lcamp = logclosure_amplitudes(vt)
    cphase = closure_phases(vt)

    uvp = VLBIFiles.load(VLBIFiles.UVData, joinpath(@__DIR__, "..", "examples", "Data", "polarized_gaussian_nogains_withdterms_withfr.uvfits"))
    dcoh = reset_mounts(VisTable(uvp, Correlations()), ARRAYFILE)

    return test_model(θ, nothing, nothing), vt, amp, lcamp, cphase, dcoh
end


function test_model(θ, _, _)
    m1 = θ.f1 * rotated(stretched(Gaussian(), θ.σ1 * θ.τ1, θ.σ1), θ.ξ1)
    m2 = θ.f2 * rotated(stretched(Gaussian(), θ.σ2 * θ.τ2, θ.σ2), θ.ξ2)
    return m1 + shifted(m2, θ.x, θ.y)
end

function test_skymodel_polarized(θ, _, metadata)
    (; lp) = metadata
    m1 = θ.f1 * rotated(stretched(Gaussian(), θ.σ1 * θ.τ1, θ.σ1), θ.ξ1)
    m2 = θ.f2 * rotated(stretched(Gaussian(), θ.σ2 * θ.τ2, θ.σ2), θ.ξ2)
    mI = m1 + shifted(m2, θ.x, θ.y)
    return PolarizedModel(mI, lp * mI, lp / 2 * mI, 0.02 * mI)
end


function test_model2(θ, _, _)
    m2 = θ.f * stretched(ExtendedRing(θ.α), θ.r * (1 + θ.τ), θ.r)
    return m2
end

function test_prior2()
    return (
        f = VLBIUniform(0.8, 1.2),
        r = VLBIUniform(μas2rad(5.0), μas2rad(30.0)),
        τ = VLBIUniform(0.0, 0.5),
        α = VLBIUniform(2.0, 8.0),
    )
end

function testimg(θ, grid, _)
    (; c, σimg) = θ
    rast = to_simplex(CenteredLR(), σimg .* c.params)
    m = ContinuousImage(rast, grid, BSplinePulse{3}())
    return m
end

function testimg_prior(grid)
    cprior = corr_image_prior(grid, μas2rad(20.0))
    prior = (c = cprior, σimg = VLBIExponential(0.1))
    return prior
end

function testimg_add(θ, grid, _)
    (; c, σimg, fg) = θ
    rast = to_simplex(CenteredLR(), σimg .* c.params)
    m = ContinuousImage((1 - fg) .* rast, grid, BSplinePulse{3}())
    g = fg * modify(Gaussian(), Stretch(μas2rad(500.0)))
    return m + g
end

function testimg_add_prior(grid)
    cprior = corr_image_prior(grid, μas2rad(20.0))
    prior = (c = cprior, σimg = VLBIExponential(0.1), fg = VLBIUniform(0.0, 1.0))
    return prior
end


function test_prior()
    return (
        f1 = VLBIUniform(0.8, 1.2),
        σ1 = VLBIUniform(μas2rad(1.0), μas2rad(40.0)),
        τ1 = VLBIUniform(0.35, 0.65),
        ξ1 = VLBIUniform(-π / 2, π / 2),
        f2 = VLBIUniform(0.3, 0.7),
        σ2 = VLBIUniform(μas2rad(1.0), μas2rad(40.0)),
        τ2 = VLBIUniform(0.35, 0.65),
        ξ2 = VLBIUniform(-π / 2, π / 2),
        x = VLBIUniform(-μas2rad(40.0), μas2rad(40.0)),
        y = VLBIUniform(-μas2rad(40.0), μas2rad(40.0)),
    )
end

# The one-channel table `t` with a second channel at frequency `ν` holding the same values.
function with_channel(t::VisTable, ν)
    a = info(t)
    ch = (; frequency = [only(a.channels.frequency), ν], width = fill(only(a.channels.width), 2))
    order = sortperm(ch.frequency)
    a2 = ArrayInfo(; stations = a.stations, channels = map(c -> c[order], ch), scans = a.scans, tracks = a.tracks, source = a.source, epoch = a.epoch)
    V, w = collect(baseimage(measurement(t))), collect(baseimage(weight(t)))
    return VisTable(t.pol, records(t), hcat(V, V), hcat(w, w); info = a2)
end

# The per-record coordinates of `t` as a records table.
function records(t::VisTable)
    d = domain(t)
    c = Comrade.ComradeBase.coords(d)
    return (; uvw = tuple.(c.u, c.v, c.w), d.Ti, d.ant1, d.ant2)
end

# `t` observed `days` days later: the same records with the epoch moved.
function later(t::VisTable, days)
    a = info(t)
    b = ArrayInfo(; a.stations, a.channels, a.scans, a.tracks, a.source, epoch = a.epoch + Dates.Day(days))
    return VisTable(t.pol, records(t), collect(baseimage(measurement(t))), collect(baseimage(weight(t))); info = b)
end
