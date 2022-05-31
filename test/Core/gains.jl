@testset "gains" begin
    _,vis, amp, lcamp, cphase = load_data()

    st = scantable(vis)
    stamp = scantable(amp)
    stcp = scantable(cphase)
    stlca = scantable(lcamp)

    gcache = GainCache(st)

    tel = stations(vis)
    gamp_prior = NamedTuple{Tuple(tel)}(ntuple(_->LogNormal(0.0, 0.1), length(tel)))
    gph_prior = NamedTuple{Tuple(tel)}(ntuple(_->Normal(0.0, 0.1), length(tel)))

    gamp = GainPrior(gamp_prior, st)
    gpha = GainPrior(gph_prior, st)

    rga = rand(gamp)
    rph = rand(gpha)
    @inferred logdensityof(gamp, rga)
    @inferred logdensityof(gpha, rph)


end
