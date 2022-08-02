using Comrade
using Distributions
using BenchmarkTools

load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
obs.add_scans()
obsavg = obs.avg_coherent(0.0, scan_avg=true)
amp = extract_amp(obsavg)
lklhd = RadioLikelihood(amp)

function model(θ)
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
post = Posterior(lklhd, prior, model)

θ = (rad= 22.0, wid= 3.0, a = 0.0, b = 0.15, f=0.8, sig = 20.0, asy=0.2, pa=π/2, x=20.0, y=20.0)
m = model(θ)

post = Posterior(lklhd, prior, model)
tpost = asflat(post)

x0 = inverse(tpost, θ)

ℓ = logdensityof(tpost)
@benchmark ℓ($x0)

using ForwardDiff
gℓ = Comrade.make_pullback(ℓ, AD.ForwardDiffBackend())
@benchmark gℓ($x0)

# Now we do the eht-imaging benchmarks

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

# This is a hack to get the objective function and its gradient
# we need to do this since the functions depend on some global variables
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
data1, sigma1, uv1, _ = ehtim.modeling.modeling_utils.chisqdata(obsavg, "amp", pol="I")
data2, sigma2, uv2, _ = ehtim.modeling.modeling_utils.chisqdata(obsavg, false, pol="I")
data3, sigma3, uv3, _ = ehtim.modeling.modeling_utils.chisqdata(obsavg, false, pol="I")

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
                "Obsdata"=>obsavg,
                "fit_pol"=>false, "fit_cpol"=>false,
                "flux"=>1.0, "alpha_flux"=>0, "fit_gains"=>false,
                "marginalize_gains"=>false, "ln_norm"=>1314.33,
                "param_init"=>pinit, "test_gradient"=>false
            )

# This is the negative log-posterior
fobj = ehtim.modeling.modeling_utils.objfunc

# This is the gradient of the negative log-posterior
gfobj = ehtim.modeling.modeling_utils.objgrad

using BenchmarkTools
@benchmark fobj($pinit)

@benchmark gfobj($pinit)
