# Introduction to the VLBI Imaging Problem

Very-long baseline interferometry (VLBI) is capable of taking the highest resolution images in the world, achieving angular resolutions of ~20 μas. In 2019, the first-ever image of a black hole was produced by the Event Horizon Telescope (EHT). However, while the EHT has unprecedented resolution, it is also a sparse interferometer. As a result, the sampling in the uv or Fourier space of the image is incomplete. This incompleteness makes the imaging problem uncertain. Namely, infinitely many images are possible, given the data. `Comrade` is a
imaging/modeling package that aims to quantify this uncertainty using
Bayesian inference.

If we denote visibilities by `V` and the image structure/model by `I`, `Comrade` will then compute the posterior or the probability of an image given the visibility data or in an equation

```math
p(I|V) = \frac{p(V|I)p(I)}{p(V)}.
```

Here ``p(V|I)`` is known as the likelihood and describes the probability distribution of the data given some image `I`. The prior ``p(I)`` encodes prior knowledge of the image structure. This prior includes distributions of model parameters and even the model itself. Finally, the denominator ``p(V)`` is a normalization term and is known as the marginal likelihood or evidence and can be used to assess how well particular models fit the data.

Therefore, we must specify the likelihood and prior to construct our posterior. Below we provide a brief description of the likelihoods and models/priors that `Comrade` uses. However, if the user wants to see how everything works first, they should check out the [Making an Image of a Black Hole](@ref) tutorial.

## Likelihood

Following TMS[^TMS], we note that the likelihood for a single complex visibility at baseline ``u_{ij}, v_{ij}`` is

```math 
p(V_{ij} | I) = (2\pi \sigma^2_{ij})^{-1/2}\exp\left(-\frac{| V_{ij} - g_ig_j^*\tilde{I}_{ij}(I)|^2}{2\sigma^2_{ij}}\right).
```

In this equation, ``\tilde{I}`` is the Fourier transform of the image ``I``, and ``g_{i,j}`` are complex numbers known as gains. The gains arise due to atmospheric and telescope effects and corrupt the incoming signal. Therefore, if a user attempts to model the complex visibilities, they must also model the complex gains. An example showing how to model gains in `Comrade` can be found in [Stokes I Simultaneous Image and Instrument Modeling](@ref).

Modeling the gains can be computationally expensive, especially if our image model is simple. For instance, in `Comrade`, we have a wide variety of geometric models. These models tend to have a small number of parameters and are simple to evaluate. Solving for gains then drastically increases the amount of time it takes to sample the posterior. As a result, part of the typical EHT analysis[^M87P6][^SgrAP4] instead uses closure products as its data. The two forms of closure products are:

  - Closure Phases,
  - Log-Closure Amplitudes.

Closure Phases ``\psi`` are constructed by selecting three baselines ``(i,j,k)`` and finding the argument of the bispectrum:

```math
    \psi_{ijk} = \arg V_{ij}V_{jk}V_{ki}.
```

Similar log-closure amplitudes are found by selecting four baselines ``(i,j,k,l)`` and forming the closure amplitudes:

```math
    A_{ijkl} = \frac{ |V_{ij}||V_{kl}|}{|V_{jk}||V_{li}|}.
```

Instead of directly fitting closure amplitudes, it turns out that the statistically better-behaved data product is the log-closure amplitude. 

The benefit of fitting closure products is that they are independent of complex gains, so we can leave them out when modeling the data. However, the downside is that they effectively put uniform improper priors on the gains[^Blackburn], meaning that we often throw away information about the telescope's performance. On the other hand, we can then view closure fitting as a very conservative estimate
about what image structures are consistent with the data. Another downside of using closure products is that their likelihoods are complex. In the high-signal-to-noise limit, however, they do reduce to Gaussian likelihoods, and this is the limit we are usually in for the EHT. For the explicit likelihood `Comrade` uses, we refer the reader to [appendix F](https://iopscience.iop.org/article/10.3847/2041-8213/ac6736#apjlac6736app6) in paper IV of the first Sgr A* EHT publications[^SgrAP4]. The computational implementation of these likelihoods can be found in [VLBILikelihoods.jl](https://github.com/ptiede/VLBILikelihoods.jl).

## Prior Model

`Comrade` has included a large number of possible models (see [Comrade API](@ref) for a list). These can be broken down into two categories:

  1. Parametric or *geometric models*
  2. Non-parametric or *image models*

`Comrade`'s geometric model interface is different from other EHT modeling packages because we don't directly provide fully formed models. Instead, we offer simple geometric models, which we call **primitives**. These primitive models can then be **modified** and combined to form complicated 
image structures. For more information, we refer the reader to [Model Interface](@ref).

Additionally, we include an interface to Bayesian imaging methods, where we directly fit a rasterized image to the data. These models are highly flexible and assume very little about the image structure. In that sense, these methods are an excellent way to explore the data first and see what kinds of image structures are consistent with observations. For an example of how to fit an image model to closure products, we refer the reader to the other tutorial included in the docs.

## References

[^TMS]: Thompson, A., Moran, J., Swenson, G. (2017). Interferometry and Synthesis in Radio Astronomy (Third). Springer Cham
[^M87P6]: Event Horizon Telescope Collaboration, (2022). First M87 Event Horizon Telescope Results. VI. The Shadow and Mass of the Central Black Hole. ApJL 875 L6 [doi](https://doi.org/10.3847/2041-8213/ab1141)
[^SgrAP4]: Event Horizon Telescope Collaboration, (2022). First Sagittarius A* Event Horizon Telscope Results. IV. Variability, Morphology, and Black Hole Mass. ApJL 930 L15 [arXiv](https://doi.org/10.3847/2041-8213/ac6736)
[^Blackburn]: Blackburn, L., et. al. (2020). Closure statistics in interferometric data. ApJ, 894(1), 31.
