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
 - name: Center for Astrophysics | Harvard & SmithsonianS
   index: 2
date: 9 Feb 2022
bibliography: paper.bib

---

# Summary

- Mention VLBI
- Fully Bayesian modeling of VLBI data (Bayesian modeling and inference of VLBI data)
- Black holes mention them and 



`Comrade` is a package for Bayesian modeling and inference using radio interferometric measurements. `Comrade` is targeted to very-long-baseline-interferometry researchers and has been widely used across the EHT. `Comrade` can model complicated source structures, either by composing simple geometric image models, or by physical modeling of emission in a curved spacetime. Through its Bayesian framework `Comrade` provides a posterior estimates for these geometric and physical parameters.
 

<!-- `Comrade` is written in Julia and is designed to take advantage of Julia's differentiable programming, and high-performance nature. Julia was chosen to allow for end-users to incorporate their own models, while maintaining high-performance.  -->

<!-- "and then to be fit directly to a flexible range of interferometric data products such as interferometric visibilities, visibility amplitudes, and closure quantities."  -->

# Statement of need

Radio interferometric measurement provide the highest resolution images every produced, culminating in the first image of a black hole. However, because interferometers such as the EHT, only provide sparse sampling in the Fourier domain, accurate quantification of uncertainty in the high-dimensional image space is both imperative and numerically demanding. 

`Comrade` is a Julia package designed to efficiently quantify the uncertainty in radio imaging, while remaining dynamic and modular. Julia is a natural language for this, due to its python-esque syntax and C-like performance. Additionally, due to Julia's extensive auto-differnetiation libraries, `Comrade` can provide gradients to optimization and sampling algorithms. 

<!-- "Comrade is designed to be especially useful for VLBI studies of black holes, with tailored model classes..." 

This ability to differentiate models natively is important when considering complicated source morphologies as expected . In such datasets models with large numbers of parameters are needed. Therefore, gradient accelerated optimization/sampling algorithms are needed to efficiently solve the problem. 

`Comrade` was designed to be used by radio astronomers during analysis of VLBI data. It has already been used in a number of analyses for the Event Horizon Telescope that will soon lead to publications. 

"It has already been used in a number of analyses for the Event Horizon Telescope that will soon lead to publications." << "It is used widely across the EHT Collaboration for analysis and interpretation." 

# Mathematics

`Comrade` was written with VLBI in mind. In particular, we assume that the image size $\ll 1$rad. In this limit the ideal visibilities are given by the Fourier transform of the image specific intensity $I_\nu(\alpha, \beta)$ :
$$
    \mathcal{V}(u,v) = \int e^{2\pi i (u\alpha + v\beta)}I(\alpha, \beta)\mathrm{d}\alpha\mathrm{d}\beta.
$$

`Comrade` provides an interface to quickly specify an image structure and it's resulting Fourier transform. The general problem of VLBI is then inverting this relation. That is, moving from a set of measured visibilities $V(u, v)$ to an image structure. This is complicated by the fact that visibility measurements are sparse. This makes the inverse problem degenerate, and variety of source structures are possible. 

To solve this problem, `Comrade` views the problem as a Bayesian inverse problem. Therefore, `Comrade` provides a variety of source model classes and likelihood functions applicable for VLBI data analysis.
Comrade itself does not explicity include any optimizers or samplers to find the optimal images. This is by design. Selecting the appropriate sampler often depends on the best the data set, image model, etc. Instead `Comrade` makes it easy to construct a log posterior density and then fit it with your preferred optimizer. 
For instance if a user wants to use nested sampling to fit the problem they can do:
 -->
```julia
using Comrade
using Distributions
using ForwardDiff
using Pathfinder
using AdvancedHMC

# load eht-imaging we use this to load eht data
load_ehtim()
obs = ehtim.obsdata.load_uvfits("data.uvfits")
# remove zero baselines
obsflg = obs.flag_uvdistance(uvmin=0.1e9)
# Extract visibility amplitudes and closure phases
dlcamp = extract_lcamp(obs; count="min")
dcphase = extract_cphase(obs, count="min")

# form the likelihood
lklhd = RadioLikelihood(damp, dcphase)

# build the model here we fit a ring with a azimuthal brightness variation and a Gaussian
function model(θ)
  (; radius, width, α, β, fgauss, σ, τ1, ξ1, x1, y1) = θ
  ring = f1*smoothed(stretched(MRing(α, β), radius, radius), width)
  g1 = (1-f1)*shifted(rotated(stretched(Gaussian(), σ1/sqrt(1-τ), σ2*sqrt(1-τ)), ξ1), x1, y1)
  return ring + g1
end

# defines my priors
prior = ( 
          radius = Uniform(μas2rad(10.0), μas2rad(30.0)),
          width = Uniform(μas2rad(1.0), μas2rad(20.0)),
          α = Uniform(-0.5, 0.5),
          β = Uniform(-0.5, 0.5),
          f1 = Uniform(0.0, 1.0),
          σ1 = Uniform(μas2rad(1.0), μas2rad(40.0)),
          τ1 = Uniform(0.0, 0.75),
          ξ1 = Uniform(-π/2, π/2),
          x1 = Uniform(-μas2rad(60.0), μas2rad(60.0))
          y1 = Uniform(-μas2rad(60.0), μas2rad(60.0))
        )
# Now form my posterior
post = Posterior(lklhd, prior, model)
# transform posterior to (-∞,∞) space 
tpost = asflat(post)
logp(x) = logdensity(tpost, x)



```


However, `Comrade` can be used by itself. Some examples of this are given in its documentation.

# Exosystem

`Comrade` is broken up into two main github repos: Comrade.jl and ComradeBase.jl. Comrade is the package non-developers can use, and includes plotting utilities, data management, and optimization/sampling interfaces. `ComradeBase` defines the inferface that Comrade uses, and is designed to be a low dependency package used mostly by developers. This allows someone to include their own radio source models without having to depend on the full `Comrade` package.

To solve the Bayesian inverse problem `Comrade` has a number of small interfaces with different Julia probabilistic programming packages. Currently the most developed version is `ComradeSoss`, which defines a minimal interface between `Comrade` and `Soss`. In addition, an interface with `Turing` and `BAT` is planned. 


# Similar Packages

- <span style="font-variant:small-caps;">Themis</span>: A C++ parameter estimation package used by the EHT. It is currently a private GitHub repo
- `eht-imaging`: Python general purpose EHT imaging package. It currently has a modeling submodule. Requires hand written gradients.
- `Galifrey`: Python modeling package that uses emcee and no gradients for model fitting
- `eht-dmc`: Python Bayesian imaging package that also fits calibration systematics by solving the RIME [@Hamaker:1996].

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements


# References