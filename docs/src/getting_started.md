# Getting Started

To get started with ROSE

```julia
using ROSE

#load ehtim
load_ehtim()

# To load some data we can use
obs = ehtim.obsdata.load_uvfits("FILENAME")

# Extracts the visibility amplitudes
dvis   = extract_vis(obs)
dcp    = extract_cphase(obs)
dlcamp = extract_lcamp(obs)

# Just want the array configuration?
ac = arrayconfig(data)

# Want to construct a elliptical Gaussian
m1 = rotated(stretched(Gaussian(), 20.0, 10.0), π/4)
m2 = ExtendedRing(20.0)
mtot = m1+m2

# Now construct a model image wrapper
mimage = modelimage(mtot)

#Now evaluate the visibilities
vis = visibility(mimage, ac)
cphase = closure_phase(mimage, )
```