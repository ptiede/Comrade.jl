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
# we need to do this since the functions depend on some global variables
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
# 298 us +/- 7.7

# This is the gradient of the negative log-posterior
gfobj = eh.modeling.modeling_utils.objgrad
# 1.3 ms


