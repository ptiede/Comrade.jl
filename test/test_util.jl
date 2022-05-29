function load_data()
    obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "test_data.uvfits"))
    obs.add_scans()
    obsavg = obs.avg_coherent(0.0, scan_avg=true)

    m = ehtim.model.Model()
    m = m.add_gauss(1.0, μas2rad(40.0), μas2rad(20.0), π/3, 0.0, 0.0)
    m = m.add_gauss(0.5, μas2rad(20.0), μas2rad(10.0), π/6, μas2rad(30.0), μas2rad(30.0))
    obsm = m.observe_same_nonoise(obsavg)

    vis = extract_vis(obsm)
    amp = extract_amp(obsm)
    lcamp = extract_lcamp(obsm)
    cphase = extract_cphase(obsm, cut_trivial=true)

    return m, vis, amp, lcamp, cphase
end


function test_model(θ)
    m1 = θ.f1*rotated(stretched(Gaussian(), θ.σ1*θ.τ1, θ.σ1), θ.ξ1)
    m2 = θ.f2*rotated(stretched(Gaussian(), θ.σ2*θ.τ2, θ.σ2), θ.ξ2)
    return m1 + shifted(m2, θ.x, θ.y)
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
