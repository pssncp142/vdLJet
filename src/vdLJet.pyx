# cython: language_level=3

cimport cython
from cython.parallel import prange

from libc.math cimport sin, cos, tan, sqrt, pow, exp, isnan, M_PI
from libc.stdio cimport printf

import numpy as np
cimport numpy as np

DTYPE = np.double

ctypedef np.double_t DTYPE_t

cdef double c = 2.997e8 # m/s

# Returns Lorentz factor
# beta_b : bulk velocity in units of c
cdef double _Lorentz_factor(double beta_b) nogil except 0:
    return sqrt(1./(1.-beta_b**2))

# Returns Doppler factor
# beta_b : bulk velocity in units of c
# inc : inclication angle in units of pi
# mode : -1 or 1 (approaching, receding)
cdef double _Doppler_factor(double beta_b, double inc, int mode) nogil except 0:
    if mode == -1:
        return 1./(_Lorentz_factor(beta_b)*(1-beta_b*cos(inc)))
    elif mode == 1:
        return 1./(_Lorentz_factor(beta_b)*(1+beta_b*cos(inc)))

# Equation 11
cdef double _beta_exp(double beta_b, double inc, double tan_theta_obs) nogil:
    cdef double Lorentz_f = _Lorentz_factor(beta_b)
    cdef double bulk_f = sqrt(pow(Lorentz_f,2)*(1-pow(beta_b*cos(inc),2))-1.)
    return tan_theta_obs*bulk_f

# Equation 4 
cdef double _p_from_tau_0(double tau_0) nogil:
    return 3./(2.*tau_0)*(exp(tau_0)-tau_0-1)

# Equation 5
cdef double _radius(double delta_t_r, double R_0, double beta_exp) nogil:
    return R_0 + beta_exp*c*delta_t_r

# Equation 12
cdef double _t_0_from_tej(double t_ej, double R_0, double beta_exp) nogil except 0:
    return t_ej+R_0/beta_exp/c

cdef double _tau_nu(double delta_t_r, double nu, double R_0, double beta_exp,
                    double tau_0, double nu_0, double p) nogil:
    cdef double radius_t = _radius(delta_t_r, R_0, beta_exp)
    cdef double radius_t_factor = pow(radius_t/R_0, -(2.*p+3.))
    cdef double nu_factor = pow(nu/nu_0, -(p+4.)/2.)
    return tau_0*radius_t_factor*nu_factor

cdef double _flux_nu_ej_rest(double t, double nu, double F_0, double R_0, double beta_exp,
                    double tau_0, double nu_0, double t_0, double p) nogil:
    cdef double delta_t_r = (t-t_0)
    cdef double radius_t = _radius(delta_t_r, R_0, beta_exp)
    cdef double radius_t_factor = pow(radius_t/R_0, 3.)
    cdef double nu_factor = pow(nu/nu_0, 5./2.)
    cdef double tau_nu_t = _tau_nu(delta_t_r, nu, R_0, beta_exp, tau_0, nu_0, p)
    cdef double tau_nu_t_factor = (1.-exp(-tau_nu_t))/(1.-exp(-tau_0))
    cdef double res = F_0*radius_t_factor*nu_factor*tau_nu_t_factor
    if isnan(res):
        res = 0
    
    return res

cpdef double _from_F_0_to_R_0(double F_0, double d, double tau_0):
    """
    Fill out
    """
    return sqrt(F_0*pow(d, 2)/M_PI/(1-exp(-tau_0)))

def flux_nu_ej(np.ndarray[DTYPE_t, ndim=1] t, double nu, dict config, int mode, bint verbose=False):
    """
    flux_nu_ej(time, nu, config, mode)

    Calculate observed jet evolution based on van der Laan model.

    Paramaters
    ----------
    - time : (N,) ndarray
        Timestamps to calculate the flux
    - nu : double
        Frequency to calculate lightcurve
    - config : dict
        adasd
    - mode : int 
        (-1, 1) approaching/receding jet
    
    Returns
    -------
    - flux : (N,) ndarray
        Flux density calculated at input timestamps
    """

    if verbose:
        print("Flux_0     : %.2e Jy" % config['F_0'])
        print("p          : %.2f   " % _p_from_tau_0(config['tau_0']))
        print("tau_0      : %.2f   " % config['tau_0'])
        print("R_0        : %.2e m " % config['R_0'])
        print("t_ej       : %.2e s " % config['t_ej'])
        print("nu_0       : %.2e Hz" % config['nu_0'])
        print("beta_b     : %.2e c " % config['beta_b'])
        print("theta_obs  : %.1f deg " % config['theta_obs'])
        print("inc        : %.1f deg " % config['inc'])

    if (mode != -1 and mode != 1):
        return None

    return flux_nu_ej_raw(t, nu, config['F_0'], config['R_0'], config['beta_b'], config['theta_obs'], 
                        config['inc'], config['tau_0'], config['t_ej'], config['nu_0'], mode)

@cython.boundscheck(False)
@cython.wraparound(False) 
cpdef flux_nu_ej_raw(np.ndarray[DTYPE_t, ndim=1] t, double nu, double F_0, double R_0, double beta_b, 
                    double theta_obs, double inc, double tau_0, double t_ej,
                    double nu_0, int mode):

    cdef int i
    cdef int sz = t.shape[0]
    cdef np.ndarray res = np.zeros([sz], dtype=DTYPE)
    cdef double[:] res_view = res
    cdef double[:] t_view = t
    cdef double tan_theta_obs = tan(theta_obs/180.*M_PI)

    cdef double tmp_beta_exp = _beta_exp(beta_b, inc/180.*M_PI, tan_theta_obs)
    cdef double doppler_f = _Doppler_factor(beta_b, inc/180.*M_PI, mode)
    cdef double p = _p_from_tau_0(tau_0)
    cdef double t_0 = _t_0_from_tej(t_ej/doppler_f, R_0, tmp_beta_exp)

    for i in prange(sz, nogil=True):
        res_view[i] = doppler_f**3*_flux_nu_ej_rest(t_view[i], nu*doppler_f, F_0, R_0, tmp_beta_exp, tau_0, nu_0, t_0, p)

    return res





    


    
