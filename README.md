# ROSE

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ptiede.github.io/ROSE.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ptiede.github.io/ROSE.jl/dev)
[![Build Status](https://github.com/ptiede/ROSE.jl/workflows/CI/badge.svg)](https://github.com/ptiede/ROSE.jl/actions)
[![Coverage](https://codecov.io/gh/ptiede/ROSE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ptiede/ROSE.jl)

Radio Observation Sampling Exercises

This is a alpha version of a potential EHT/ngEHT modeling/imaging software written in pure Julia. The objectives are in order of importance

1. Simple to define and adjust models and parameterizations
2. Everything must be differentiable with Julia's autodifferentation libraries
3. Speed. The code must be fast, with proper SIMD, GPU, and distributed computing support

# Installation
To install ROSE.jl you must use Julia's package manager. For example you can launch your Julia session with then type `]` to move into Pkg mode. Once in pkg mode type
```julia
add "https://github.com/ptiede/ROSE.jl"
```
to install the main branch of ROSE.jl. Eventually this will be push to the official Julia repo's but the name of the package may change before that happens.


# Roadmap

1. Settle on model interface
2. Design a better data maniuplation interface and native uvfits reading (or spin up the EHTIM.jl package)
4. Better method to store and save array and image configurations to speed up everything
5. Full gain and polarization implemented, including DTERM modeling ala DMC, and approximate methods such as pseudo-marginal, particle MCMC methods, or laplace approximations
   - Include a trait the defines whether the model is image and visibility analytic, and have this pass through composite models 
5. Better sampling and model update interface 
   - Probably will use Accessors.jl and some other tool to do this
6. Write model adjoints when needed for additional speed (and define a may to define AD backends like Zygote, Enzyme, Diffractor, etc.)
7. Hook up into PPL's?
   - Already have something like this for Soss.jl at [ROSESoss.jl](https://github.com/ptiede/ROSESoss.jl)
   - Add something for Turing and Gen?
8. GPU support when its a good idea (maybe for gain solves?)

