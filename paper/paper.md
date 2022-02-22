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

To show the simplicity of `Comrade`'s interface we will reproduce results from @EHTCVI in under 50 lines of code which finishes in under 2 minutes.
![](code.pdf)

![Image of M 87 from `Comrade`](blackhole.png)


# Similar Packages

- [<span style="font-variant:small-caps;">Themis</span>](themis): A C++ parameter estimation package used by the EHT. It is currently a private GitHub repo. 
- `eht-imaging`: Python general purpose EHT imaging package. It currently has a modeling submodule.
- `Galifrey`: Python modeling package that uses emcee as it's backbone.
- `eht-dmc`: Python Bayesian imaging package that also fits calibration systematics by solving the RIME [@Hamaker:1996].

# Acknowledgements

PT thanks Michael Johnson, Dom Pesce, Lindy Blackburn, and Avery Broderick for their helpful discussions related to the development of this package and VLBI in general.

# References