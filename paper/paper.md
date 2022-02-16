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

`Comrade` is a package for modeling different radio astronomy source structures. It has been designed to allow for complicated source structures to be constructed from simple geometric models, and then be compared to data. `Comrade` is targeted to very-long-baseline-interferometry researchers. It will be valuable for analyzing what source structures are supported by the data, in a Bayesian modeling framework. By modeling the VLBI imaging problem as a Bayesian inverse problem, `Comrade` can provide uncertainty quantification of image structures, that are not typically possible with typical tools. 

`Comrade` is written in Julia and is designed to take advantage of Julia's differentiable programming, and 

Very-long-baseline interfeometry (VLBI) is capable of producing the highest resolution images ever produced. In 2019 the first ever image of a black hole was presented in 

# Statement of need

`Comrade` is a Julia package for modeling radio emission. Julia is a high performance dynamic language that bridges the "two-language problem". This allows for users to easily extend `Comrade` while maintaining the speed of a lower-level language (e.g. C). The API for `Comrade` was designed to allow the user to easily add their own emission models whether the image has an analytic Fourier transform. Additionally, `Comrade` is a differentiable emission modeling library. This allows for the use of gradient information to be passed to optimization and Bayesian inference algorithms.

This ability to differentiate models natively is important when considering complicated source morphologies commonly seen in VLBI datasets. In such datasets models with large numbers of parameters are needed. Therefore, gradient accelerated optimization/sampling algorithms are needed to efficiently solve the problem. 

`Comrade` was designed to be used by radio astronomers during analysis of VLBI data. It has already been used in a number of analyses for the Event Horizon Telescope that will soon lead to publications. 

# Mathematics

`Comrade` was written with VLBI in mind. In particular, we assume that the image size $\ll 1$rad. In this limit the ideal visibilities are given by the Fourier transform of the image specific intensity $I_\nu(\alpha, \beta)$ :
$$
    \mathcal{V}(u,v) = \int e^{2\pi i (u\alpha + v\beta)}I(\alpha, \beta)\mathrm{d}\alpha\mathrm{d}\beta.
$$

`Comrade` provides an interface to quickly specify an image structure and it's resulting Fourier transform. The general problem of VLBI is then inverting this relation. That is, moving from a set of measured visibilities $V(u, v)$ to an image structure. This is complicated by the fact that visibility measurements are sparse. This makes the inverse problem degenerate, and variety of source structures are possible. 

To solve this problem, `Comrade` uses views the problem as a Bayesian inverse problem. Therefore, `Comrade` provides a variety of source model classes and likelihood functions applicable for VLBI data analysis.

To solve the Bayesian inverse problem `Comrade` has a number of small interfaces with different Julia probabilistic programming packages. Currently the most developed version is `ComradeSoss`, which defines a minimal interface between `Comrade` and `Soss`. In addition, an interface with `Turing` and `BAT` is planned. 

However, `Comrade` can be used by itself. Some examples of this are given in its documentation.




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