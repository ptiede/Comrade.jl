---
title: 'Comrade: Composable Modeling of Radio Emission'
tags:
  - julia
  - astronomy
  - radio astronomy
  - vlbi
  - black holes
authors:
  - name: Paul Tiede
    orcid: 0000-0003-3826-5648
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Black Hole Initiative at Harvard University
   index: 1
 - name: Center for Astrophysics | Harvard & Smithsonian
   index: 2
date: 9 Feb 2022
bibliography: paper.bib

---

# Summary

`Comrade` is a Bayesian modeling package, targeted for very-long-baseline interferometry (VLBI) and written in the Julia[^1] programming language [@bezanson2015julia]. `Comrade` aims at producing VLBI image of black holes and active galactic nuclei. Furthermore, it focuses on providing uncertainty quantification of image and physical source properties, such as a black hole's accretion state. The package has already been widely used within the Event Horizon Telescope Collaboration and will be useful for expert and novice VLBI researchers.

[^1]: https://julialang.org


# Statement of need

Radio interferometric measurements provide the highest resolution images ever produced, culminating in the first image of a black hole [@EHTCI; @EHTCIV; @EHTCVI]. However, producing VLBI images is not straightforward.
An ideal VLBI array samples the Fourier transform of the image, $I$:

$$
V(u,v) = \int I(\alpha, \beta) e^{-2\pi i (u\alpha + v\beta)}d\alpha d\beta.
$$

In general, VLBI data sets provide an incomplete sampling of $V(u_i, v_i)$ in the Fourier domain. Therefore, VLBI images are inherently uncertain and quantifying this uncertainty is fundamental to the VLBI imaging problem. This quantification is especially significant for the Event Horizon Telescope, which typically has only 5-8 distinct observing sites. To model this uncertainty, `Comrade` uses Bayesian inference and casts VLBI imaging as a Bayesian inverse problem.

Given the diverse nature of AGN and black holes, `Comrade` includes geometric models, such as Gaussians, disks, rings, crescents to extract relevant features. For non-parametric modeling/imaging, `Comrade` includes a rasterized image model similar to the one described in @themaging. Finally,  `Comrade`'s flexible model interface, enables direct physical modeling of an accretion flow in curved spacetime.

Bayesian inference is numerically demanding relative to traditional VLBI imaging and modeling methods. Traditionally these computational demands have required writing large sections of the code in a lower-level language, e.g., `C`/`C++`. However, this approach comes at a productivity cost to the end-user and makes it difficult for researchers to add their own models. `Comrade` solves this problem by using the Julia programming language. Julia was designed to solve this two-language problem by having `C`-like performance while maintaining a Python-esque syntax and programming experience [@juliafast]. Julia achieves these features using a just-ahead-of-time compiler, and code specialization based on multiple dispatch. The effectiveness of this approach has been demonstrated in, e.g, the Celeste project [@celeste], where Julia was the first dynamically typed language to break the petaflops barrier.

As Julia is a differentiable programming language, most `Comrade` models are natively differentiable. Utilizing gradient information to explore the parameter space quickly will be necessary to find reasonable image structures as image complexity grows. For images of the central black hole of AGN, this is the norm due to complicated accretion structure. `Comrade` is therefore well equipped to deal with the Bayesian VLBI imaging problem, even for large VLBI arrays.

To sample from the posterior `Comrade` has interfaces to nested sampling algorithms, `NestedSamplers` [@ns] and `AdvancedHMC` [@AHMC] by default. Moreover, users can specify their own samplers. To make it easy to use `Comrade` with other posterior samplers, the `Comrade` includes functionality that transforms the posterior density from the parameter space to $\mathbb{R}^n$, as well as the unit hypercube. This functionality is needed for Hamiltonian Monte Carlo and nested sampling respectively.

As a result of `Comrade`'s design, it makes it easy to quickly produce posteriors of VLBI data. Here we show an example program that reproduces results from @EHTCVI. Namely and produces posterior of the image structure of the black hole in M 87^[On an Intel i5-7200U chip, this finishes in 165s].

```julia
using Comrade
using Distributions
using ComradeAHMC
using ComradeOptimization
using OptimizationBBO
using Plots
# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
file = "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"
obs = ehtim.obsdata.load_uvfits(file)
obs.add_scans()
# kill 0-baselines since we don't care about 
# large scale structure and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
# form the likelihood
lklhd = RadioLikelihood(dlcamp, dcphase)
# build the model: here we fit a ring with a azimuthal 
# brightness variation and a Gaussian
function model(params)
  (;rad, wid, a, b, f, sig, asy, pa, x, y) = params
  ring = f*smoothed(stretched(MRing((a,), (b,)), rad, rad), wid)
  g = (1-f)*shifted(rotated(stretched(Gaussian(), sig*asy, sig), pa), x, y)
  return ring + g
end
# define the priors
uas2rad = pi/180.0/3600/1e6
prior = (
          rad = Uniform(uas2rad*(10.0), uas2rad*(30.0)),
          wid = Uniform(uas2rad*(1.0), uas2rad*(10.0)),
          a = Uniform(-0.5, 0.5), b = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          sig = Uniform(uas2rad*(1.0), uas2rad*(40.0)),
          asy = Uniform(0.0, 0.75),
          pa = Uniform(0.0, 1pi),
          x = Uniform(-uas2rad*(80.0), uas2rad*(80.0)),
          y = Uniform(-uas2rad*(80.0), uas2rad*(80.0))
        )
# Now form the posterior
post = Posterior(lklhd, prior, model)
# We will use HMC to sample the posterior.
# First we will find a reasonable starting location using Optimization
# to have nice bounds we first transform to the unit hypercube
tpost = ascube(post)
ndim = dimension(tpost)
f = OptimizationFunction(tpost)
prob = OptimizationProblem(
            f, rand(ndim), nothing;
            lb=fill(1e-2, ndim), ub = fill(0.99, ndim)
            )
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000)
# transform the solution back to regular space
xopt = transform(tpost, sol.u)
# Comrade is all about uncertainty quantification so now let's find the posterior!
# To do this we will use the `AdvancedHMC` package or rather its interface to Comrade.
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric), 3000; nadapts=2000, init_params=xopt)
```

![Output of the above code. The image is a random posterior draw for an image of M 87.](blackhole.png)


# Similar Packages

- `eht-imaging` [@chael2018]: Python general purpose EHT imaging package. It currently has a modeling submodule.
- `eht-dmc` [@dmc]: Python Bayesian polarized imaging package that also fits calibration systematics by solving the radio interferometry measurement equation [@Hamaker].
- `Galifray` [@gal]: Python modeling package that uses emcee as its sampler.
- `InterferometricModels` [@in]: Recent Julia radio astronomy package with some similar features to `Comrade`.
- `NIFTy`[@nifty]: Python-based bayesian imaging only package that uses Gaussian processes to model images including VLBI images. 
- <span style="font-variant:small-caps;">Themis</span> [@themis]: A C++ parameter estimation package used by the EHT. It is currently a private repository for the EHT.


# Acknowledgements

P.T. thanks Michael Johnson, Dom Pesce, Lindy Blackburn, and Avery Broderick for their helpful discussions related to the development of this package and VLBI in general.

This work was supported by the Black Hole Initiative at Harvard University, which is funded by grants from the John Templeton Foundation and the Gordon and Betty
Moore Foundation to Harvard University. Additional support was provided by Perimeter Institute for Theoretical Physics. Research at Perimeter Institute is supported by the Government of Canada through the Department of Innovation, Science
and Economic Development Canada and by the Province
of Ontario through the Ministry of Economic Develop-
ment, Job Creation and Trade.  PT's research is also supported by the National Science Foundation (AST-1935980 and AST-2034306) and additional financial support from
the Natural Sciences and Engineering Research Council of
Canada through a Discovery Grant.

# References