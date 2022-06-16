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
