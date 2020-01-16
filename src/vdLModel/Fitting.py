import vdLJet
import numpy as np

def flux_nu_ej_rec(time, nu, i, params):
    v = params.valuesdict()

    return vdLJet.flux_nu_ej_raw(time, nu, v['F_0_%d' % (i+1)], v['R_0_%d' % (i+1)], v['beta_b_%d' % (i+1)],
                                           v['theta_obs_%d' % (i+1)], v['inc_%d' % (i+1)],
                                           v['tau_0_%d' % (i+1)], v['t_ej_%d' % (i+1)], v['n_0_%d' % (i+1)], 1)

def flux_nu_ej_app(time, nu, i, params):
    v = params.valuesdict()

    return vdLJet.flux_nu_ej_raw(time, nu, v['F_0_%d' % (i+1)], v['R_0_%d' % (i+1)], v['beta_b_%d' % (i+1)],
                                           v['theta_obs_%d' % (i+1)], v['inc_%d' % (i+1)],
                                           v['tau_0_%d' % (i+1)], v['t_ej_%d' % (i+1)], v['n_0_%d' % (i+1)], -1)

def flux_nu_ej_both_one(time, nu, i, params):

    return flux_nu_ej_rec(time, nu, i, params) + flux_nu_ej_app(time, nu, i, params)

def flux_nu_ej_both_all(time, nu, n, params):

    res = np.zeros(time.shape)
    for i in range(n):
        res += flux_nu_ej_both_one(time, nu, i, params)

    return res

def flux_nu_ej_both_all_resid(params, time, data, err, nu, n):

    model = flux_nu_ej_both_all(time, nu, n, params)
    return (data-model)/err

def make_params(n, values, mins, maxs, varys):
    pass
 