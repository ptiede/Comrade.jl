using Comrade
using Distributions
using BenchmarkTools
using Pyehtim

# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "..", "examples", "Data", "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obsavg = scan_average(obs)
amp = extract_table(obsavg, VisibilityAmplitudes())

function model(θ, p)
    (;rad, wid, a, b, f, sig, asy, pa, x, y) = θ
    ring = f*smoothed(stretched(MRing((a,), (b,)), μas2rad(rad), μas2rad(rad)), μas2rad(wid))
    g = (1-f)*shifted(rotated(stretched(Gaussian(), μas2rad(sig)*asy, μas2rad(sig)), pa), μas2rad(x), μas2rad(y))
    return ring + g
end
prior = (
          rad = Uniform(10.0, 30.0),
          wid = Uniform(1.0, 10.0),
          a = Uniform(-0.5, 0.5), b = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          sig = Uniform((1.0), (60.0)),
          asy = Uniform(0.0, 0.9),
          pa = Uniform(0.0, 1π),
          x = Uniform(-(80.0), (80.0)),
          y = Uniform(-(80.0), (80.0))
        )
# Now form the posterior
skym = SkyModel(model, prior, imagepixels(μas2rad(150.0), μas2rad(150.0), 128, 128))

θ = (rad= 22.0, wid= 3.0, a = 0.0, b = 0.15, f=0.8, sig = 20.0, asy=0.2, pa=π/2, x=20.0, y=20.0)
m = model(θ, nothing)

post = VLBIPosterior(skym, amp)
tpost = asflat(post)

x0 = prior_sample(tpost)

ℓ = logdensityof(tpost)
@benchmark ℓ($x0)
# 38.1 μs
using Enzyme
@benchmark Enzyme.gradient(Enzyme.Reverse, $(Const(tpost)), $x0)
#107.3 μs

# Now we do the eht-imaging benchmarks
