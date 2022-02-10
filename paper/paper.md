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

Very-long-baseline interfeometry (VLBI) is capable of producing the highest resolution images ever produced.  

# Statement of need

`Comrade` is a Julia package for modeling radio emission. Julia is a high performance dynamic language that bridges the "two-language problem". This allows for users to easily extend `Comrade` while maintaining the speed of a lower-level language (e.g. C). The API for `Comrade` was designed to allow the user to easily add their own emission models whether the image has an analytic Fourier transform. Additionally, `Comrade` is a differentiable emission modeling library. This allows for the use of gradient information to be passed to optimization and Bayesian inference algorithms. 

`Comrade` was designed to be used by radio astronomers during analysis of VLBI data. It has already been used in a number of analyses for the Event Horizon Telescope that soon lead to publications. Additionally, it has been used for 

# Mathematics

Currently `Comrade` assumes that the field of view of small enough that ignoring telescope systematics, are given by the Fourier transform
$$
    \mathcal{V}(u,v) = \int e^{2\pi i (u\alpha + v\beta)}I(\alpha, \beta)\mathrm{d}\alpha\mathrm{d}\beta.
$$



`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics


# Similar Packages

- <span style="font-variant:small-caps;">Themis</span>: A C++ parameter estimation package used by the EHT. It is currently a private GitHub repo
- `eht-imaging`: Python general purpose EHT imaging package. It currently has a modeling submodule. Requires hand written gradients.
- `Galifrey`: Python modeling package that uses emcee and not gradients for model fitting
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

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References