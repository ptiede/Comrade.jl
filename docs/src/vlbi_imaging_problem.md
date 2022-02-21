# Getting Started

## Introduction to the VLBI Imaging Problem

Very-long baseline interferometry (VLBI) is capable of taking the highest resolution images in the world, achieving angular resolutions of ~20 μas. In 2019 the first ever image of a black hole was produced by the Event Horizon Telescope (EHT). However, the EHT is a very sparse interferometer. This makes imaging uncertain. This is because VLBI doesn't measure an actual image, but rather the Fourier transform of the on-sky image. Futhermore, the EHT only have ~7 dishes at any one time. This makes the Fourier sampling very sparse/incomplete. As a result, there is no single image that can explain the data. Traditional imaging methods get around this with curated imaging algorithms, such as DIFMAP and regularized-maximum-likelihood. However, these methods still ignore the inherent uncertainty in the image reconstruction. `Comrade` instead views the imaging problem as a Bayesian inverse problem. 



Namely, there isn't a single image capable of explaining the source structure but a family. `Comrade` is a Julia software package that views the imaging problem as a classical Bayesian inverse problem.

Due to the high resolution VLBI data can be quite complicated. As such there is a variety of possible source structures. As such, `Comrade` provides a number of 

To get started with Comrade

```julia
using Comrade

#load ehtim
load_ehtim()

# To load some data we can use
obs = ehtim.obsdata.load_uvfits("FILENAME")

# Extract various data products from the EHTIM object
dvis   = extract_vis(obs)
dcp    = extract_cphase(obs)
dlcamp = extract_lcamp(obs)

# Just want the array configuration?
ac = arrayconfig(data)

# Want to construct a elliptical Gaussian
m1 = rotated(stretched(Gaussian(), 20.0, 10.0), π/4)
m2 = ExtendedRing(20.0, 10.0)
mtot = m1+m2

# plot the model
plot(mtot)

# Lets make an image!
img = intensitymap(mtot, 80.0, 80.0, 128, 128)

# plot the image
plot(img)

# Now construct a model image wrapper
# Since Extended ring doesn't have a simple FT we need to specify an image to hold it.
# In the future this will be automated.
mimage = modelimage(mtot, img)


#Now evaluate the visibilities
vis = visibilities(mimage, ac)
```
