# Getting Started

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
m2 = ExtendedRing(20.0)
mtot = m1+m2

# Lets make an image!
img = intensitymap(mtot, 80.0, 80.0, 128, 128)

# Now construct a model image wrapper
# Since Extended ring doesn't have a simple FT we need to specify an image to hold it.
# In the future this will be automated.
mimage = modelimage(mtot, img)


#Now evaluate the visibilities
vis = visibility(mimage, ac)
```
