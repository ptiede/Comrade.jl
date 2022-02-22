---
title: 'Comrade: Composable Modeling of Radio Emission'
tags:
  - Julia
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

`Comrade` is a Bayesian modeling package targeted for very-long-baseline interferometry written in the Julia programming language [@bezanson2015julia]. `Comrade` aims at providing VLBI image reconstructions of black holes and active galactic nuclei. Furthermore, it focuses on providing uncertainty quantification of image and physical source properties, such as the accretion state of the black hole. The package has already been widely used within the Event Horizon Telescope Collaboration and will generally be useful for VLBI researchers and students.


# Statement of need

Radio interferometric measurements provide the highest resolution images every produced, culminating in the first image of a black hole [@EHTCI; @EHTCIV; @EHTCVI]. However, producing images is not straightforward.
An ideal VLBI telescope measures the Fourier transform of the image, $I$

$$
V(u,v) = \int I(\alpha, \beta) e^{2\pi i (u\alpha + v\beta)}d\alpha d\beta.
$$

In general, VLBI data sets provide an incomplete sampling of $u_i, v_j$ in the Fourier domain. Therefore, it is impossible to construct a single image $I$ given the measurements $V(u_i, v_j)$. Instead there is infinitely many possible images that are "consistent" with the data. This is especially relevant for the Event Horizon Telescope, typically has only 4-5 unique baselines before aperture synthesis. Therefore, to quantify what image structures are robust, a notion of uncertainty is required. To model this uncertainty, `Comrade` uses Bayesian inference and views the VLBI imaging problem as a classical Bayesian inverse problem. However, Bayesian inference is numerically demanding relative to traditional VLBI imaging methods. This motivates the need for a  high-performance language. However, having a dynamic and interactive language is also important for day-to-day use and user friendliness. Previous packages had separated this into two languages, where a high-performance language (e.g., C, C++) was used for Bayesian inference, and then all inspection and visualizations were written in a separate high-level language (e.g., python). `Comrade` solves this problem by using the Julia programming language. Julia was designed to solve this two-language problem and it able to achieve `C`-like speeds, while maintaining a python-esque syntax and programming experience.

As Julia is a differentiable programming language, all `Comrade` models are natively differentiable. This is unique for an EHT modeling library where either gradients have to be hand-coded, or are calculated using finite difference. The use of gradient information is imperative for VLBI modeling. Namely, as data sets get larger so will the number of model parameters needed to describe the data. Utilizing gradient information to quick explore the parameter space will then be necessary to find reasonable image structures. `Comrade` is therefore well equipped to deal with the VLBI imaging problem in the next few years.

To sample from the posterior `Comrade` has some interfaces to nested sampling algorithms `NestedSamplers.jl` and `AdvancedHMC` [@AHMC] by default, but it is easy for users to use different samplers. To make it easy to use `Comrade` with other samplers, the package includes functionality that transforms the posterior density from the parameter space $P$ to $\mathbb{R}^n$, as well as the unit hypercube, which is needed for Hamiltonian Monte Carlo and nested sampling respectively.

As a result of `Comrade`'s design, it makes it easy to quickly produces posteriors of VLBI data. Here we show an example program that reproduces some of the results from [EHTCVI], and produces posterior of the image structure of the black hole in M 87

```julia
using Comrade
using Distributions
using Pathfinder
using AdvancedHMC
using Plots
# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
file = "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"
obs = ehtim.obsdata.load_uvfits(file)
obs.add_scans()
# kill 0-baselines since we don't care about 
# large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs; count="min")
dcphase = extract_cphase(obs, count="min")
# form the likelihood
lklhd = RadioLikelihood(dlcamp, dcphase)
# build the model here we fit a ring with a azimuthal 
#brightness variation and a Gaussian
function model(θ)
  (;rad, wid, a, b, f, sig, asy, pa, x, y) = θ
  ring = f*smoothed(stretched(MRing((a,), (b,)), rad, rad), wid)
  g = (1-f)*shifted(rotated(stretched(Gaussian(), sig*asy, sig), pa), x, y)
  return ring + g
end
# define the priors
uas2rad = π/180.0/3600/1e6
prior = (
          rad = Uniform(uas2rad*(10.0), uas2rad*(30.0)),
          wid = Uniform(uas2rad*(1.0), uas2rad*(10.0)),
          a = Uniform(-0.5, 0.5), b = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          sig = Uniform(uas2rad*(1.0), uas2rad*(40.0)),
          asy = Uniform(0.0, 0.75),
          pa = Uniform(0.0, 1π),
          x = Uniform(-uas2rad*(80.0), uas2rad*(80.0)),
          y = Uniform(-uas2rad*(80.0), uas2rad*(80.0))
        )
# Now form the posterior
post = Posterior(lklhd, prior, model)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
q, phi, _ = multipathfinder(post, 100)
# now we sample using hmc
metric = DiagEuclideanMetric(dimension(post))
chain, stats = sample(post, HMC(;metric), 2000; 
                      nadapts=1000, init_params=phi[1])
# plot a draw from the posterior
plot(model(chain[end]))    
```

![Image of M 87 from `Comrade`](blackhole.png)


# Similar Packages

- `eht-imaging` [@chael2018]: Python general purpose EHT imaging package. It currently has a modeling submodule.
- `eht-dmc` [@dmc]: Python Bayesian imaging package that also fits calibration systematics by solving the RIME [@Hamaker].
- `Galifray` [@gal]: Python modeling package that uses emcee as it's backbone.
- `InterferometricModels` [@in]: Recent Julia radio astronomy package with some similar features to `Comrade`
- [<span style="font-variant:small-caps;">Themis</span>](@themis): A C++ parameter estimation package used by the EHT. It is currently a private GitHub repo.


# Acknowledgements

PT thanks Michael Johnson, Dom Pesce, Lindy Blackburn, and Avery Broderick for their helpful discussions related to the development of this package and VLBI in general.

# References