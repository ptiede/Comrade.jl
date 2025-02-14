# Benchmarks

`Stoked` was partially designed with performance in mind. Solving imaging inverse problems is traditionally very computationally expensive, especially since Stoked uses Bayesian inference. To benchmark `Stoked` we will compare it to two of the most common modeling or imaging packages within the EHT:

- [eht-imaging](https://github.com/achael/eht-imaging/)
- [Themis](https://iopscience.iop.org/article/10.3847/1538-4357/ab91a4)

`eht-imaging`[^1] or `ehtim` is a Python package that is widely used within the EHT for its imaging and modeling interfaces. It is easy to use and is commonly used in the EHT. However, to specify the model, the user must specify how to calculate the model's complex visibilities **and** its gradients, allowing eht-imaging's modeling package to achieve acceptable speeds.

Themis is a C++ package focused on providing Bayesian estimates of the image structure. In fact, `Stoked` took some design cues from `Themis`. Themis has been used in various EHT publications and is the standard Bayesian modeling tool used in the EHT. However, `Themis` is quite challenging to use and requires a high level of knowledge from its users, requiring them to understand makefile, C++, and the MPI standard. Additionally, Themis provides no infrastructure to compute gradients, instead relying on finite differencing, which scales poorly for large numbers of model parameters. 

## Benchmarking Problem

For our benchmarking problem, we analyze a situation very similar to the one explained in  Namely, we will consider fitting 2017 M87 April 6 data using an m-ring and a single Gaussian component. Please see the end of this page to see the code we used for `Stoked` and `eht-imaging`.

## Results

All tests were run using the following system

```julia
Julia Version 1.10.3
Python Version 3.10.12
Stoked Version 0.10.0
eht-imaging Version 1.2.7
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
```

Our benchmark results are the following:

| | Stoked (micro sec) | eht-imaging (micro sec) | Themis (micro sec)|
|---|---|---|---|
| posterior eval (min) | 31.1  | 445  | 55  |
| posterior eval (mean) | 31.8  | 476  | 60  |
| grad posterior eval (min) |  104 (Enzyme) | 1898  | 1809  |
| grad posterior eval (mean) |  107 (Enzyme) | 1971 |  1866  |

Therefore, for this test we found that `Stoked` was the fastest method in all tests. For the posterior evaluation we found that Stoked is > 10x faster than `eht-imaging`, and 2x faster then `Themis`. For gradient evaluations we have `Stoked` is > 15x faster than both `eht-imaging` and `Themis`.

[^1]: Chael A, et al. *Interferometric Imaging Directly with Closure Phases* 2018 ApJ 857 1 arXiv:1803/07088

## Code

### Julia Code

```julia
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "..", "examples", "Data", "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obsavg = scan_average(obs)
amp = extract_table(obsavg, VisibilityAmplitudes())

function model(θ, p)
    (;rad, wid, a, b, f, sig, asy, pa, x, y) = θ
    ring = f*smoothed(stretched(MRing((a,), (b,)), μas2rad(rad), μas2rad(rad)), μas2rad(wid))
    g = (1-f)*shifted(rotated(stretched(Gaussian(), μas2rad(sig)*asy, μas2rad(sig)), pa), μas2rad(x), μas2rad(y))
    return ring + g
end
prior = (
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
# Now form the posterior
skym = SkyModel(model, prior, imagepixels(μas2rad(150.0), μas2rad(150.0), 128, 128))

θ = (rad= 22.0, wid= 3.0, a = 0.0, b = 0.15, f=0.8, sig = 20.0, asy=0.2, pa=π/2, x=20.0, y=20.0)
m = model(θ, nothing)

post = VLBIPosterior(skym, amp)
tpost = asflat(post)

x0 = prior_sample(tpost)

ℓ = logdensityof(tpost)
@benchmark ℓ($x0)
# 31.1 μs
using Enzyme
@benchmark Enzyme.gradient(Enzyme.Reverse, $(Const(ℓ)), $x0)
# 104 μs
```

### eht-imaging Code

```python
import ehtim as eh
import numpy as np
import os
obs = eh.obsdata.load_uvfits(os.path.join("examples", "Data", "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
obsavg = obs.avg_coherent(0.0, scan_avg=True)

meh = eh.model.Model()
meh = meh.add_thick_mring(F0=0.8,
                    d=2*22.0*eh.RADPERUAS,
                    alpha=2*np.sqrt(2*np.log(2))*eh.RADPERUAS*3.0,
                    x0 = 0.0,
                    y0 = 0.0,
                    beta_list=[0.0+0.15j]
                    )
meh = meh.add_gauss(F0=1-0.8,
                    FWHM_maj=2*np.sqrt(2*np.log(2))*eh.RADPERUAS*(3.0),
                    FWHM_min=2*np.sqrt(2*np.log(2))*eh.RADPERUAS*(3.0)*0.2,
                    PA = np.pi/2,
                    x0 = eh.RADPERUAS*(20.0),
                    y0 = eh.RADPERUAS*(20.0)
                    )

preh = meh.default_prior()
preh[0]["F0"] = {"prior_type": "flat", "min" : 0.0, "max" : 1.0}
preh[0]["d"] = {"prior_type": "flat", "min" : eh.RADPERUAS*(20.0), "max" : eh.RADPERUAS*(60.0)}
preh[0]["alpha"] = {"prior_type": "flat", "min" : eh.RADPERUAS*(2.0), "max" : eh.RADPERUAS*(25.0)}
preh[0]["x0"] = {"prior_type": "fixed"}
preh[0]["y0"] = {"prior_type": "fixed"}

preh[1]["F0"] = {"prior_type": "flat", "min" : 0.0, "max" : 1.0}
preh[1]["FWHM_maj"] = {"prior_type": "flat", "min" : eh.RADPERUAS*(2.0), "max" : eh.RADPERUAS*(120.0)}
preh[1]["FWHM_min"] = {"prior_type": "flat", "min" : eh.RADPERUAS*(2.0), "max" : eh.RADPERUAS*(120.0)}
preh[1]["x0"] = {"prior_type": "flat", "min" : -eh.RADPERUAS*(40.0), "max" : eh.RADPERUAS*(40.0)}
preh[1]["y0"] = {"prior_type": "flat", "min" : -eh.RADPERUAS*(40.0), "max" : eh.RADPERUAS*(40.0)}
preh[1]["PA"] = {"prior_type": "flat", "min" : -np.pi, "max" : np.pi}

# This is a hack to get the objective function and its gradient
# we need to do this since the functions depend on some global ehtim variables
transform_param = eh.modeling.modeling_utils.transform_param
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

# make the python param map and use optimize so we flatten the parameter space.
pmap, pmask = eh.modeling.modeling_utils.make_param_map(meh, preh, "scipy.optimize.dual_annealing", fit_model=True)
trial_model = meh.copy()

# get initial parameters
n_params, pinit = make_paraminit(pmap, meh, trial_model, preh)

# make data products for the globdict
data1, sigma1, uv1, _ = eh.modeling.modeling_utils.chisqdata(obsavg, "amp", pol="I")
data2, sigma2, uv2, _ = eh.modeling.modeling_utils.chisqdata(obsavg, True, pol="I")
data3, sigma3, uv3, _ = eh.modeling.modeling_utils.chisqdata(obsavg, True, pol="I")
data4, sigma4, uv4, _ = eh.modeling.modeling_utils.chisqdata(obsavg, True, pol="I")

# now set the ehtim modeling globdict

eh.modeling.modeling_utils.globdict = {"trial_model" : trial_model,
                "d1" : "amp", "d2" : False, "d3" : False, "d4" : False,
                "pol1" : "I", "pol2" : "I", "pol3" : "I", "pol4" : "I",
                "data1" : data1, "sigma1" : sigma1, "uv1" : uv1, "jonesdict1" : None,
                "data2" : data2, "sigma2" : sigma2, "uv2" : uv2, "jonesdict2" : None,
                "data3" : data3, "sigma3" : sigma3, "uv3" : uv3, "jonesdict3" : None,
                "data4" : data3, "sigma4" : sigma3, "uv4" : uv3, "jonesdict4" : None,
                "alpha_d1" : 0, "alpha_d2" : 0, "alpha_d3" : 0, "alpha_d4" : 0,
                "n_params" :  n_params, "n_gains" : 0, "n_leakage" : 0,
                "model_prior" : preh, "param_map" : pmap, "param_mask" : pmask,
                "gain_prior" : None, "gain_list" : [], "gain_init" : None,
                "fit_leakage" : False, "leakage_init" : [], "leakage_fit" : [],
                "station_leakages" : None, "leakage_prior" : None,
                "show_updates" : False, "update_interval" : 1,
                "gains_t1" : None, "gains_t2" : None,
                "minimizer_func" : "scipy.optimize.dual_annealing",
                "Obsdata" : obsavg,
                "fit_pol" : False, "fit_cpol" : False,
                "flux" : 1.0, "alpha_flux" : 0, "fit_gains" : False,
                "marginalize_gains" : False, "ln_norm" : 1314.33,
                "param_init" : pinit, "test_gradient" : False
}

# This is the negative log-posterior
fobj = eh.modeling.modeling_utils.objfunc
%timeit fobj(pinit)
# 298 us +/- 7.7

# This is the gradient of the negative log-posterior
gfobj = eh.modeling.modeling_utils.objgrad
%timeit gfobj(pinit)
# 1.3 ms
```
