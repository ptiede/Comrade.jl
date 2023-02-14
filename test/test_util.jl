function load_data()
    obs = load_ehtim_uvfits(joinpath(@__DIR__, "test_data.uvfits"))
    obs.add_scans()
    obsavg = scan_average(obs)

    m = ehtim.model.Model()
    m = m.add_gauss(1.0, μas2rad(40.0), μas2rad(20.0), π/3, 0.0, 0.0)
    m = m.add_gauss(0.5, μas2rad(20.0), μas2rad(10.0), π/6, μas2rad(30.0), μas2rad(30.0))
    obsm = m.observe_same_nonoise(obsavg)

    vis = extract_vis(obsm)
    amp = extract_amp(obsm)
    lcamp = extract_lcamp(obsm)
    cphase = extract_cphase(obsm, cut_trivial=true)
    dcoh = extract_coherency(obsm)

    return m, vis, amp, lcamp, cphase, dcoh
end


function test_model(θ)
    m1 = θ.f1*rotated(stretched(Gaussian(), θ.σ1*θ.τ1, θ.σ1), θ.ξ1)
    m2 = θ.f2*rotated(stretched(Gaussian(), θ.σ2*θ.τ2, θ.σ2), θ.ξ2)
    return m1 + shifted(m2, θ.x, θ.y)
end

function test_model2(θ, metadata)
    (; alg, fovx, fovy, nx, ny) = metadata
    m2 = θ.f*stretched(ExtendedRing(θ.α), θ.r*(1+θ.τ), θ.r)
    return modelimage(m2; fovx, fovy, nx, ny, alg, thread=true)
end

function test_prior2()
    return (f = Uniform(0.8, 1.2),
            r = Uniform(μas2rad(5.0), μas2rad(30.0)),
            τ = Uniform(0.0, 0.5),
            α = Uniform(2.0, 8.0)
            )
end

function test_prior()
    return (f1=Uniform(0.8, 1.2),
             σ1 = Uniform(μas2rad(1.0), μas2rad(40.0)),
             τ1 = Uniform(0.35, 0.65),
             ξ1 = Uniform(-π/2, π/2),
             f2 = Uniform(0.3, 0.7),
             σ2 = Uniform(μas2rad(1.0), μas2rad(40.0)),
             τ2 = Uniform(0.35, 0.65),
             ξ2 = Uniform(-π/2, π/2),
             x = Uniform(-μas2rad(40.0), μas2rad(40.0)),
             y = Uniform(-μas2rad(40.0), μas2rad(40.0))
            )
end
