using Pkg; Pkg.activate(@__DIR__)
using Comrade
using Distributions
using ComradeOptimization
using ComradeDynesty
using ComradeAHMC
using OptimizationBBO
using OptimizationOptimJL
using Plots
using StatsBase

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_101_lo_hops_netcal_StokesI.uvfits"))
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = scan_average(obs.add_fractional_noise(0.01))
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs; cut_trivial=true)
# form the likelihood
# build the model here we fit a ring with a azimuthal
#brightness variation and a Gaussian

function asymgaussian(f, σ, τ, ξ, x, y)
    return f*shifted(rotated(stretched(Gaussian(), σ, σ*τ), ξ), x, y)
end


function model(θ)
    (;rout, ψ, τ, ξ, σc, slash, fl,
      f1, σ1, τ1, ξ1, x1, y1,
      f2, σ2, τ2, ξ2, x2, y2,
      f3, σ3, τ3, ξ3, x3, y3,
      fG, σG
       ) = θ

    rinner = (1-ψ)*rout
    shift = τ*(rout-rinner)

    cres = rotated(ConcordanceCrescent(rout, rinner, shift, slash), ξ)
    floor = rotated(shifted(stretched(Disk(), rinner, rinner), shift, zero(typeof(shift))), ξ)

    g1 = asymgaussian(f2, σ1, τ1, ξ1, x1, y1)
    g2 = asymgaussian((1-f2)*f3, σ2, τ2, ξ2, x2, y2)
    g3 = asymgaussian(1-(f2 + (1-f2)*f3), σ3, τ3, ξ3, x3, y3)
    compact = (1-f1)*(smoothed(cres*(1-fl) + fl*floor, σc)) + f1*(g1 + g2 + g3)
    gL = stretched(Gaussian(), σG, σG)
    return (1-fG)*compact + fG*gL
end
# define the priors
prior = (
          rout = Uniform(μas2rad(10.0), μas2rad(30.0)),
          ψ = Uniform(0.0, 1.0),
          τ = Uniform(0.0, 1.0),
          ξ = Uniform(0.0, 2π),
          σc = Uniform(μas2rad(0.1), μas2rad(15.0)),
          slash = Uniform(0.0, 1.0),
          fl = Uniform(0.0, 1.0),
          fG = Uniform(0.0, 1.0),
          σG = Uniform(μas2rad(500.0), μas2rad(2000.0)),
          f1 = Uniform(0.0, 1.0),
          σ1 = Uniform(μas2rad(1.0), μas2rad(50.0)),
          τ1 = Uniform(0.1, 1.0),
          ξ1 = Uniform(0.0, 1π),
          x1 = Uniform(-μas2rad(100.0), μas2rad(100.0)),
          y1 = Uniform(-μas2rad(100.0), μas2rad(100.0)),
          f2 = Uniform(0.0, 1.0),
          σ2 = Uniform(μas2rad(1.0), μas2rad(50.0)),
          τ2 = Uniform(0.1, 1.0),
          ξ2 = Uniform(0.0, 1π),
          x2 = Uniform(-μas2rad(100.0), μas2rad(100.0)),
          y2 = Uniform(-μas2rad(100.0), μas2rad(100.0)),
          f3 = Uniform(0.0, 1.0),
          σ3 = Uniform(μas2rad(1.0), μas2rad(50.0)),
          τ3 = Uniform(0.1, 1.0),
          ξ3 = Uniform(0.0, 1π),
          x3 = Uniform(-μas2rad(100.0), μas2rad(100.0)),
          y3 = Uniform(-μas2rad(100.0), μas2rad(100.0)),
        )

# Now form the posterior
lklhd = RadioLikelihood(model, dlcamp, dcphase)
post = Posterior(lklhd, prior)
# We will use HMC to sample the posterior.
# First we will find a reasonable starting location using GalacticOptim
# For optimization we need to specify what transform to use. Here we will transform to
# the unit hypercube
tpost = ascube(post)
ndim = dimension(tpost)
f = OptimizationFunction(tpost, Optimization.AutoForwardDiff())
sols = map(1:1) do _
    prob = OptimizationProblem(f, rand(ndim), nothing, lb=fill(0.01, ndim), ub = fill(0.99, ndim))
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=150_000, show_trace=true)
end

lmin, imin = findmin(getproperty.(sols, :minimum))
sol = sols[imin]
# Now let's get the Laplace approximation since it is cheap!
#prob = OptimizationProblem(f, sol.u, nothing)
#ldist = laplace(prob, LBFGS(); show_trace=true, maxiters=5_000, g_tol=1e-2)

# transform the solution back to regular space
xopt = transform(tpost, sol)
c2 = chi2(model(xopt), dlcamp, dcphase)/(length(dlcamp) + length(dcphase) - dimension(post))
# Let's see if the best fit looks reasonable by plotting normalized residuals for the
# log-closure amplitudes
residual(model(xopt), dlcamp)
residual(model(xopt), dcphase)


# we can also plot the best fit model or maximum likelihood estimate (MLE)
X,Y = imagepixels(μas2rad(200.0), μas2rad(200.0), 256, 256, 0.0, 0.0)
img = intensitymap(model(xopt), (;X, Y))
plot(img, title="MLE M87")

using ComradeDynesty
a1 = NestedSampler(;ndim=dimension(post), nlive=1000, sample="rwalk")
a2 = DynamicNestedSampler(dimension(post))

chain, stats = sample(post, a1; dlogz=0.01)


# chain, stats = sample(post, a2)

echain = sample(chain, Weights(stats.weights), 1000) |> Table


const PRIOR = (
    rad = Uniform(10.0, 30.0),# Radius of primary component
    inner_rad_ratio = Uniform(0.1, 0.99),# Radius ratio of Crescent cutout
    τm = Uniform(0.1, 1.),# Ratio of Stretched MRing major axis with respect to minor axis
    ξm = Uniform(0, 2π),# Angular position of Stretched MRing major axis or Crescent PA
    rg = Uniform(2.0, 40.0),# Radial size of Gaussians
    τg = Uniform(0.1, 0.99),# Ratio of Elliptical Gaussian major axis with respect to minor axis
    ξg = Uniform(-π/2, π/2),# Angular position of Elliptical Gaussian major axis
    xg = Uniform(-50.0, 50.0),# Gaussian horizontal displacement with respect to primary
    yg = Uniform(-50.0, 50.0),# Gaussian vertical displacement with respect to primary
    f = Uniform(0., 1.),# Flux ratio between primary and secondary components
    fg = Uniform(0.51, 1),# Flux ratio between Gaussians (Forces ordering based on flux)
    shift_ratio = Uniform(0.01, 0.99),# Crescent cutout offset from center
    slash= Uniform(0.0, 1.0),# Crescent slash
    floor = Uniform(0.0, 1.0),# Crescent Emission floor flux ratio
    width = Uniform(0.1, 50.0),# Gaussian blur width
    width_ratio =  Uniform(0.1, 0.99),# Double MRing Gaussian blur ratio between primary and secondary MRings. (Used in 2 MRing models)
    amp = Uniform(0.0, 0.5),# Amplitude of MRing modes
    phase = Uniform(0.0, 2π),# Phase of MRing modes
)

function Two_Elliptical_Gaussians(prior)
    params = (rg=2, ξg=2, τg=2, xg=1, yg=1, fg=1)
    gauss1, _ = Elliptical_Gaussian(prior)
    gauss2, _ = Elliptical_Gaussian(prior)

    two_ellip_gauss(θ) = begin
        (;rg, ξg, τg, xg, yg, fg) = θ
        prm1 = (rg=rg[1], ξg=ξg[1], τg=τg[1])
        prm2 = (rg=rg[2], ξg=ξg[2], τg=τg[2])

        fg*gauss1(prm1) + (1-fg)*(shiftedμas(gauss2(prm2), xg, yg))
    end

    return (model=two_ellip_gauss, prior=_model_prior(params, prior))
end

function BlurredAndSlashed_Crescent(prior)
    params = (rad=1, inner_rad_ratio=1, shift_ratio=1, ξm=1, width=1, slash=1)

    bl_sl_crescent(θ) = begin
        (;rad, inner_rad_ratio, shift_ratio, ξm, width, slash) = θ

        radμas = μas2rad(rad)
        shift = shift_ratio*(radμas *(1-inner_rad_ratio))
        inner_rad = inner_rad_ratio*radμas

        @chain ConcordanceCrescent(radμas, inner_rad ,shift, slash) begin
            smoothedμas(_, width)
            rotated(_, ξm)
        end
    end

    return (model=bl_sl_crescent, prior=_model_prior(params, prior))
end

function BlurredAndSlashed_Crescent_And_2_Elliptical_Gaussians(prior)

    params = (rad=1, inner_rad_ratio=1, shift_ratio=1, ξm=1, ξg=2, width=1, slash=1,rg=2, τg=2, xg=2, yg=2, f=1,fg=1)

    bl_sl_crescent_2_ell_gauss(θ) = begin
        (;rg, τg, ξg, xg, yg, f, fg) = θ
        b_s_cresc, _ = BlurredAndSlashed_Crescent(prior)
        el_gauss, _ = Two_Elliptical_Gaussians(prior)

        gauss_prm = (rg=rg, τg=τg, ξg=ξg, xg=xg[2], yg=yg[2], fg=fg)
        comp1 = f*b_s_cresc(θ)
        comp2 = (1-f)*shiftedμas(el_gauss(gauss_prm), xg[1], yg[1])

        comp1 + comp2
    end

    return (model=bl_sl_crescent_2_ell_gauss, prior=_model_prior(params, prior))
end


function BlurredAndSlashed_Crescent_And_Emission_Floor_And_3_Elliptical_Gaussians(prior)
    params = (rad=1, inner_rad_ratio=1, shift_ratio=1, ξm=1, ξg=2, width=1, slash=1, rg=2, τg=2, floor=1, xg=2, yg=2, f=1, fg=1)

    bl_sl_crescent_floored_2_ell_gauss(θ) = begin
        (;rad, inner_rad_ratio, shift_ratio, ξm, width, floor) = θ
        b_s_cresc_2_el_gauss, _ = BlurredAndSlashed_Crescent_And_3_Elliptical_Gaussians(prior)
        inner_rad = inner_rad_ratio*rad
        shift = shift_ratio*(rad*(1-inner_rad_ratio))

        comp1 = b_s_cresc_2_el_gauss(θ)

        cresc_floor = @chain Disk() begin
                stretchedμas(_, inner_rad, inner_rad)
                shiftedμas(_, shift, 0.)
                rotated(_,ξm)
                smoothedμas(_,width[1])
        end

        comp1 + floor*cresc_floor
    end

    return (model=bl_sl_crescent_floored_2_ell_gauss, prior=_model_prior(params, prior))
end
