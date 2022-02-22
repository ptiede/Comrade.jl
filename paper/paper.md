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
    orcid: 0000-0000-0000-0000
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

- Mention VLBI
- Fully Bayesian modeling of VLBI data (Bayesian modeling and inference of VLBI data)
- Black holes mention them and 

`Comrade` is a Bayesian modeling package targeted for very-long-baseline interferometry written in the Julia programming language [@bezanson2015julia]. It's target audience is VLBI researchers and students who are interested in quantifying uncertainty in image reconstructions of VLBI data.
This is pertinent for sparse interferometer, such as the Event Horizon Telescope, where `Comrade` has been widely used.

# Statement of need

Radio interferometric measurements provide the highest resolution images every produced, culminating in the first image of a black hole [@EHTCI; @EHTCIV; @EHTCVI].
An ideal VLBI telescope measure the Fourier transform of the image, $I$

$$
V(u,v) = \int I(\alpha, \beta) e^{2\pi i (u\alpha + v\beta)}d\alpha d\beta.
$$

However, because VLBI provide an incomplete sampling of $u_i, v_j$ in the Fourier domain, it is impossible to construct a single image $I$ given the measurements $V(u_i, v_j)$. Instead there is infinitely many possible images. This is especially relevant for VLBI telescope such the Event Horizon Telescope. Imaging methods that capture this uncertainty is therefore imperative. To model this uncertainty `Comrade` uses Bayesian inference, and views the VLBI imaging problem as a classical Bayesian inverse problem. However, Bayesian inference is typically numerically demanding relative to traditional VLBI imaging methods. This motivates the use of a high-performance language. However, having a dynamic and interactive language is also important for day-to-day use and user friendliness. Previous packages had separated this into two languages, where a high-performance language was used for Bayesian inference, and then all inspection and visualizations were written in a separate high-level language like python. `Comrade` solves this problem by using the Julia programming language. Julia was partially designed to solve this two-language problem. Namely, it is possible to write high-performance code, while maintaining a python-esque syntax and programming experience.

Additionally, since Julia is a differentiable programming language, all `Comrade` models are natively differentiable. This is unique for an EHT modeling library where either gradients have to be hand-coded, or are calculated using finite difference. The use of gradient information is imperative for VLBI modeling where as data sets get larger so will observed source morphology. This implies that the number of parameters one needs to fit will grow as well. In high-dimensional settings the use of gradients is necessary to efficiently sample and explore the parameter space.

`Comrade` itself does not implement any optimization or sampling techniques. Instead it creates an interface that allows for an easy construction of the un-normalized posterior density. Additionally, it handles typically error-prone and monotonous tasks of transforming the posterior density from parameter space $P$ to $\mathbb{R}^n$ and the unit hypercube, which is needed for Hamiltonian Monte Carlo and nested sampling respectively. Additionally, `Comrade` has an interface to the probabilistic programming language `Soss` [@Soss] for further Bayesian inference automation. To sample from the posterior `Comrade` has some interfaces to nested sampling algorithms `NestedSamplers.jl` [@NS] and `AdvancedHMC` [@AHMC] by default, but it is easy for users to use different samplers.


Thanks to `Comrade`'s design and Julia it is possible to fit VLBI data with a compact user-friendly experience. Below is an example that reproduces an image of a black hole from [@EHTCVI] in under 50 lines of code and finishes in under 2 min:
```julia
using Comrade
using Distributions
using Pathfinder
using AdvancedHMC
using Plots
# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits("SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits")
obs.add_scans()
# kill 0-baselines since we don't care about 
# large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true)
# grab data products we want to fit: 
# log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs; count="min")
dcphase = extract_cphase(obs, count="min")
# form the likelihood
lklhd = RadioLikelihood(dlcamp, dcphase)
# build the model here we fit a ring with a azimuthal 
#brightness variation and a Gaussian
function model(θ)
  (;radius, width, α, β, f, σG, τG, ξG, xG, yG) = θ
  ring = f*smoothed(stretched(MRing((α,), (β,)), radius, radius), width)
  g = (1-f)*shifted(rotated(stretched(Gaussian(), σG, σG*(1+τG)), ξG), xG, yG)
  return ring + g
end
# define the priors
prior = (
          radius = Uniform(μas2rad(10.0), μas2rad(30.0)),
          width = Uniform(μas2rad(1.0), μas2rad(10.0)),
          α = Uniform(-0.5, 0.5),
          β = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          σG = Uniform(μas2rad(1.0), μas2rad(40.0)),
          τG = Uniform(0.0, 0.75),
          ξG = Uniform(0.0, 1π),
          xG = Uniform(-μas2rad(80.0), μas2rad(80.0)),
          yG = Uniform(-μas2rad(80.0), μas2rad(80.0))
        )
# Now form my posterior
post = Posterior(lklhd, prior, model)
# We will use HMC to sample the posterior, first to reduce burn in we use pathfinder
# to get a good starting location that is approximately drawn from the posterior
q, ϕ, _ = multipathfinder(post, 100)
# now we sample using hmc
ndim = dimension(post)
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(post, HMC(;metric), 2000; nadapts=1000, init_params=ϕ[1])
# plot a draw from the posterior
plot(model(chain[end]))
```

![Image of M 87 from `Comrade`](blackhole.png)


# Similar Packages

- [<span style="font-variant:small-caps;">Themis</span>](themis): A C++ parameter estimation package used by the EHT. It is currently a private GitHub repo. 
- `eht-imaging`: Python general purpose EHT imaging package. It currently has a modeling submodule.
- `Galifrey`: Python modeling package that uses emcee as it's backbone.
- `eht-dmc`: Python Bayesian imaging package that also fits calibration systematics by solving the RIME [@Hamaker:1996].

# Acknowledgements

PT thanks Michael Johnson, Dom Pesce, Lindy Blackburn, and Avery Broderick for their helpful discussions related to the development of this package and VLBI in general.

# References