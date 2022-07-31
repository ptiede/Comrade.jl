# Introduction to the VLBI Imaging Problem

Very-long baseline interferometry (VLBI) is capable of taking the highest resolution images in the world, achieving angular resolutions of ~20 μas. In 2019, the first ever image of a black hole was produced by the Event Horizon Telescope (EHT). However, while the EHT has unprecedented resolution it is also a very sparse interferometer. As a result, the sampling in the uv or Fourier space of the image is incomplete. This makes the imaging problem uncertain. Namely, infinitely many images are possible given the data. `Comrade` is a
imaging/modeling package that aims to quantify this uncertainty using
Bayesian inference.

If we denote visibilities by `V` and the image structure/model by `I`, `Comrade` will then compute the posterior or the probability of an image given the visibility data, or in an equation

```math
p(I|V) = \frac{p(V|I)p(I)}{p(V)}.
```

Here ``p(V|I)`` is known as the likelihood and describes the probability distribution of the data given some image `I`. The prior ``p(I)`` encodes prior knowledge of the image structure. This prior includes distributions of model parameters and even the model itself. Finally, the denominator ``p(V)`` is a normalization term and is known as the marginal likelihood or evidence and can be used to assess how well particular models fit the data.

The goal of `Comrade` is to calculate (or approximate) the posterior of a image model. To see how to do that please see the [Making an Image of a Black Hole](@ref) tutorial.


## The Likelihood

`Comrade` uses the traditional VLBI expression for the data likelihood. Following standard VLBI data analysis[^TMS] each complex visibility ``V_{ij}`` follows a complex normal distribution

```math 
p(V_{ij} | I) = (2\pi \sigma^2_{ij})^{-1/2}\exp\left(-\frac{| V_{ij} - \hat{V}_{ij}|^2}{2\sigma^2_{ij}}\right),
```
where ``\hat{V}`` denotes the model visibility that we are comparing to the data. This likelihood forms the basis for all data analysis. The modeled $\hat{V}$ depend specifically on the model chosen. As a result `Comrade` has implemented a large number of models to describe observations. For a list of potential models please see [Comrade API](@ref).

For other related data products, e.g. visibility amplitudes, closure products, we use the high-signal-to-noise Gaussian approximations of the distributions. For the specific expressions please see EHTC 2022 Paper IV[^SgrP4].

[^TMS]: Thompson, A., Moran, J., Swenson, G. (2017). Interferometry and Synthesis in Radio Astronomy (Third). Springer Cham
[^SgrP4]: Event Horizon Telescope Collaboration, (2022). First Sagittarius A* Event Horizon Telscope Results. IV. Variability, Morphology, and Black Hole Mass. ApJL 930 L15 [doi](10.3847/2041-8213/ac667410.3847/2041-8213/ac667210.3847/2041-8213/ac6736)

