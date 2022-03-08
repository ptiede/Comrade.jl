# Comrade

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ptiede.github.io/Comrade.jl/dev/)
[![Build Status](https://github.com/ptiede/Comrade.jl/workflows/CI/badge.svg)](https://github.com/ptiede/Comrade.jl/actions)
[![Coverage](https://codecov.io/gh/ptiede/Comrade.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ptiede/Comrade.jl)

Composable Modeling of Radio Emission

This is a alpha version of a potential EHT/ngEHT modeling/imaging software written in Julia. The objectives are in order of importance

1. Simple to define and adjust models and parameterizations
2. Everything must be differentiable with Julia's autodifferentation libraries
3. Speed. The code must be fast, with proper SIMD, GPU, and distributed computing support

# Installation
To install Comrade.jl you must use Julia's package manager. For example you can launch your Julia session with then type `]` to move into Pkg mode. Once in pkg mode type
```julia
add Comrade
```

# Roadmap

1. Settle on model interface
2. Design a better data maniuplation interface and native uvfits reading (or spin up the EHTIM.jl package)
2. ~Full gain~ and polarization implemented, including DTERM modeling ala DMC, and approximate methods such as pseudo-marginal, particle MCMC methods, or laplace approximations
3. Include a trait the defines whether the model is image and visibility analytic, and have this pass through composite models 
4. Write model adjoints when needed for additional speed (and define a may to define AD backends like Zygote, Enzyme, Diffractor, etc.)
5. Hook up into PPL's?
   - Already have something like this for Soss.jl at [ComradeSoss.jl](https://github.com/ptiede/ComradeSoss.jl)
   - Add something for Turing and Gen?
6. GPU support when its a good idea (maybe for gain solves?)

