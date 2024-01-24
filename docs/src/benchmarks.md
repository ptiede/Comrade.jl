# Benchmarks

`Comrade` was partially designed with performance in mind. Solving imaging inverse problems is traditionally very computationally expensive, especially since Comrade uses Bayesian inference. To benchmark `Comrade` we will compare it to two of the most common modeling or imaging packages within the EHT:

- [eht-imaging](https://github.com/achael/eht-imaging/)
- [Themis](https://iopscience.iop.org/article/10.3847/1538-4357/ab91a4)

`eht-imaging`[^1] or `ehtim` is a Python package that is widely used within the EHT for its imaging and modeling interfaces. It is easy to use and is commonly used in the EHT. However, to specify the model, the user must specify how to calculate the model's complex visibilities **and** its gradients, allowing eht-imaging's modeling package to achieve acceptable speeds.

Themis is a C++ package focused on providing Bayesian estimates of the image structure. In fact, `Comrade` took some design cues from `Themis`. Themis has been used in various EHT publications and is the standard Bayesian modeling tool used in the EHT. However, `Themis` is quite challenging to use and requires a high level of knowledge from its users, requiring them to understand makefile, C++, and the MPI standard. Additionally, Themis provides no infrastructure to compute gradients, instead relying on finite differencing, which scales poorly for large numbers of model parameters. 

## Benchmarking Problem

For our benchmarking problem, we analyze a situation very similar to the one explained in  Namely, we will consider fitting 2017 M87 April 6 data using an m-ring and a single Gaussian component. Please see the end of this page to see the code we used for `Comrade` and `eht-imaging`.

## Results

All tests were run using the following system

```julia
Julia Version 1.7.3
Python Version 3.10.5
Comrade Version 0.4.0
eht-imaging Version 1.2.4
Commit 742b9abb4d (2022-05-06 12:58 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
```

Our benchmark results are the following:

| | Comrade (micro sec) | eht-imaging (micro sec) | Themis (micro sec)|
|---|---|---|---|
| posterior eval (min) | 31  | 445  | 55  |
| posterior eval (mean) | 36  | 476  | 60  |
| grad posterior eval (min) |  105 (ForwardDiff) | 1898  | 1809  |
| grad posterior eval (mean) |  119 (ForwardDiff) | 1971 |  1866  |

Therefore, for this test we found that `Comrade` was the fastest method in all tests. For the posterior evaluation we found that Comrade is > 10x faster than `eht-imaging`, and 2x faster then `Themis`. For gradient evaluations we have `Comrade` is > 15x faster than both `eht-imaging` and `Themis`.

[^1]: Chael A, et al. *Inteferometric Imaging Directly with Closure Phases* 2018 ApJ 857 1 arXiv:1803/07088

## Code

### Julia Code

```julia
using Pyehtim
using Comrade
using Distributions
using BenchmarkTools
using ForwardDiff
using VLBIImagePriors
using Zygote

# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "assets/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
obs = scan_average(obs)
amp = extract_table(obs, VisibilityAmplitudes())

function model(θ)
    (;rad, wid, a, b, f, sig, asy, pa, x, y) = θ
    ring = f*smoothed(modify(MRing((a,), (b,)), Stretch(μas2rad(rad))), μas2rad(wid))
    g = modify(Gaussian(), Stretch(μas2rad(sig)*asy, μas2rad(sig)), Rotate(pa), Shift(μas2rad(x), μas2rad(y)), Renormalize(1-f))
    return ring + g
end

lklhd = RadioLikelihood(model, amp)
prior = NamedDist(
          rad = Uniform(10.0, 30.0),
          wid = Uniform(1.0, 10.0),
          a = Uniform(-0.5, 0.5), b = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          sig = Uniform((1.0), (60.0)),
          asy = Uniform(0.0, 0.9),
          pa = Uniform(0.0, 1π),
          x = Uniform(-(80.0), (80.0)),
          y = Uniform(-(80.0), (80.0))
        )

θ = (rad= 22.0, wid= 3.0, a = 0.0, b = 0.15, f=0.8, sig = 20.0, asy=0.2, pa=π/2, x=20.0, y=20.0)
m = model(θ)

post = Posterior(lklhd, prior)
tpost = asflat(post)

# Transform to the unconstrained space
x0 = inverse(tpost, θ)

# Lets benchmark the posterior evaluation
ℓ = logdensityof(tpost)
@benchmark ℓ($x0)

using LogDensityProblemsAD
# Now we benchmark the gradient
gℓ = ADgradient(Val(:Zygote), tpost)
@benchmark LogDensityProblemsAD.logdensity_and_gradient($gℓ, $x0)
```

### eht-imaging Code

```julia
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "assets/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
obs = scan_average(obs)



meh = ehtim.model.Model()
meh = meh.add_thick_mring(F0=θ.f,
                    d=2*μas2rad(θ.rad),
                    alpha=2*sqrt(2*log(2))*μas2rad(θ.wid),
                    x0 = 0.0,
                    y0 = 0.0,
                    beta_list=[0.0+θ.b]
                    )
meh = meh.add_gauss(F0=1-θ.f,
                    FWHM_maj=2*sqrt(2*log(2))*μas2rad(θ.sig),
                    FWHM_min=2*sqrt(2*log(2))*μas2rad(θ.sig)*θ.asy,
                    PA = θ.pa,
                    x0 = μas2rad(20.0),
                    y0 = μas2rad(20.0)
                    )

preh = meh.default_prior()
preh[1]["F0"] = Dict("prior_type"=>"flat", "min"=>0.0, "max"=>1.0)
preh[1]["d"] = Dict("prior_type"=>"flat", "min"=>μas2rad(20.0), "max"=>μas2rad(60.0))
preh[1]["alpha"] = Dict("prior_type"=>"flat", "min"=>μas2rad(2.0), "max"=>μas2rad(25.0))
preh[1]["x0"] = Dict("prior_type"=>"fixed")
preh[1]["y0"] = Dict("prior_type"=>"fixed")

preh[2]["F0"] = Dict("prior_type"=>"flat", "min"=>0.0, "max"=>1.0)
preh[2]["FWHM_maj"] = Dict("prior_type"=>"flat", "min"=>μas2rad(2.0), "max"=>μas2rad(120.0))
preh[2]["FWHM_min"] = Dict("prior_type"=>"flat", "min"=>μas2rad(2.0), "max"=>μas2rad(120.0))
preh[2]["x0"] = Dict("prior_type"=>"flat", "min"=>-μas2rad(40.0), "max"=>μas2rad(40.0))
preh[2]["y0"] = Dict("prior_type"=>"flat", "min"=>-μas2rad(40.0), "max"=>μas2rad(40.0))
preh[2]["PA"] = Dict("prior_type"=>"flat", "min"=>-1π, "max"=>1π)

using PyCall
py"""
import ehtim
import numpy as np
transform_param = ehtim.modeling.modeling_utils.transform_param
def make_paraminit(param_map, meh, trial_model, model_prior):
    model_init = meh.copy()
    param_init = []
    for j in range(len(param_map)):
        pm = param_map[j]
        if param_map[j][1] in trial_model.params[param_map[j][0]].keys():
            param_init.append(transform_param(model_init.params[pm[0]][pm[1]]/pm[2], model_prior[pm[0]][pm[1]],inverse=False))
        else: # In this case, the parameter is a list of complex numbers, so the real/imaginary or abs/arg components need to be assigned
            if param_map[j][1].find('cpol') != -1:
                param_type = 'beta_list_cpol'
                idx = int(param_map[j][1].split('_')[0][8:])
            elif param_map[j][1].find('pol') != -1:
                param_type = 'beta_list_pol'
                idx = int(param_map[j][1].split('_')[0][7:]) + (len(trial_model.params[param_map[j][0]][param_type])-1)//2
            elif param_map[j][1].find('beta') != -1:
                param_type = 'beta_list'
                idx = int(param_map[j][1].split('_')[0][4:]) - 1
            else:
                raise Exception('Unsure how to interpret ' + param_map[j][1])

            curval = model_init.params[param_map[j][0]][param_type][idx]
            if '_' not in param_map[j][1]:
                param_init.append(transform_param(np.real( model_init.params[pm[0]][param_type][idx]/pm[2]), model_prior[pm[0]][pm[1]],inverse=False))
            elif   param_map[j][1][-2:] == 're':
                param_init.append(transform_param(np.real( model_init.params[pm[0]][param_type][idx]/pm[2]), model_prior[pm[0]][pm[1]],inverse=False))
            elif param_map[j][1][-2:] == 'im':
                param_init.append(transform_param(np.imag( model_init.params[pm[0]][param_type][idx]/pm[2]), model_prior[pm[0]][pm[1]],inverse=False))
            elif param_map[j][1][-3:] == 'abs':
                param_init.append(transform_param(np.abs(  model_init.params[pm[0]][param_type][idx]/pm[2]), model_prior[pm[0]][pm[1]],inverse=False))
            elif param_map[j][1][-3:] == 'arg':
                param_init.append(transform_param(np.angle(model_init.params[pm[0]][param_type][idx])/pm[2], model_prior[pm[0]][pm[1]],inverse=False))
            else:
                if not quiet: print('Parameter ' + param_map[j][1] + ' not understood!')
    n_params = len(param_init)
    return n_params, param_init
"""

# make the python param map and use optimize so we flatten the parameter space.
pmap, pmask = ehtim.modeling.modeling_utils.make_param_map(meh, preh, "scipy.optimize.dual_annealing", fit_model=true)
trial_model = meh.copy()

# get initial parameters
n_params, pinit = py"make_paraminit"(pmap, meh, trial_model, preh)

# make data products for the globdict
data1, sigma1, uv1, _ = ehtim.modeling.modeling_utils.chisqdata(obs, "amp")
data2, sigma2, uv2, _ = ehtim.modeling.modeling_utils.chisqdata(obs, false)
data3, sigma3, uv3, _ = ehtim.modeling.modeling_utils.chisqdata(obs, false)

# now set the ehtim modeling globdict

ehtim.modeling.modeling_utils.globdict = Dict("trial_model"=>trial_model,
                "d1"=>"amp", "d2"=>false, "d3"=>false,
                "pol1"=>"I", "pol2"=>"I", "pol3"=>"I",
                "data1"=>data1, "sigma1"=>sigma1, "uv1"=>uv1, "jonesdict1"=>nothing,
                "data2"=>data2, "sigma2"=>sigma2, "uv2"=>uv2, "jonesdict2"=>nothing,
                "data3"=>data3, "sigma3"=>sigma3, "uv3"=>uv3, "jonesdict3"=>nothing,
                "alpha_d1"=>0, "alpha_d2"=>0, "alpha_d3"=>0,
                "n_params"=> n_params, "n_gains"=>0, "n_leakage"=>0,
                "model_prior"=>preh, "param_map"=>pmap, "param_mask"=>pmask,
                "gain_prior"=>nothing, "gain_list"=>[], "gain_init"=>nothing,
                "fit_leakage"=>false, "leakage_init"=>[], "leakage_fit"=>[],
                "station_leakages"=>nothing, "leakage_prior"=>nothing,
                "show_updates"=>false, "update_interval"=>1,
                "gains_t1"=>nothing, "gains_t2"=>nothing,
                "minimizer_func"=>"scipy.optimize.dual_annealing",
                "Obsdata"=>obs,
                "fit_pol"=>false, "fit_cpol"=>false,
                "flux"=>1.0, "alpha_flux"=>0, "fit_gains"=>false,
                "marginalize_gains"=>false, "ln_norm"=>1314.33,
                "param_init"=>pinit, "test_gradient"=>false
            )

# This is the negative log-posterior
fobj = ehtim.modeling.modeling_utils.objfunc

# This is the gradient of the negative log-posterior
gfobj = ehtim.modeling.modeling_utils.objgrad

# Lets benchmark the posterior evaluation
@benchmark fobj($pinit)

# Now we benchmark the gradient
@benchmark gfobj($pinit)
```
