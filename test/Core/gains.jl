using Distributions

@testset "gains" begin
    _,vis, amp, lcamp, cphase = load_data()

    st = scantable(vis)
    stamp = scantable(amp)
    stcp = scantable(cphase)
    stlca = scantable(lcamp)

    gcache = GainCache(st)


    θ = (f1 = 1.0,
         σ1 = 1.2e-10,
         τ1 = 0.25,
         ξ1 = -0.42,
         f2 = 0.5, σ2 = 8.0e-12,
         τ2 = 0.4,
         ξ2 = 0.0,
         x = 4.0e-11,
         y = 1.5e-10)

    m = test_model(θ)

    tel = stations(vis)
    gamp_prior = NamedTuple{Tuple(tel)}(ntuple(_->LogNormal(0.0, 0.1), length(tel)))
    gph_prior = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))

    gamp = GainPrior(gamp_prior, st)
    gpha = GainPrior(gph_prior, st)

    ga = fill(1.0, size(rand(gamp))...)
    gm = GainModel(gcache, ga, m)

    ac = arrayconfig(vis)
    cpac = arrayconfig(cphase)
    lcac = arrayconfig(lcamp)
    @test visibilities(gm, ac) ≈ visibilities(m, ac)
    @test amplitudes(gm, ac) ≈ amplitudes(m, ac)
    @test closure_phases(gm, cpac) ≈ closure_phases(m, cpac)
    @test logclosure_amplitudes(gm, lcac) ≈ logclosure_amplitudes(m, lcac)


    rga = rand(gamp)
    rph = rand(gpha)
    @inferred logdensityof(gamp, rga)
    @inferred logdensityof(gpha, rph)
end
