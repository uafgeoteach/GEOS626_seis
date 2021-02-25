import numpy as np
from scipy.integrate import solve_ivp
from stress_disp_tor import stress_disp_tor
import spshell_config

def surf_stress(f):
    """ Python adaptation of surf_stress.m by Carl Tape
        Coding by Amanda McPherson, Dec 2020
        
        Integrates using the RK45 method, and calls stress_disp_tor to calculate the
        derivatives of displacement and stress
        
        INPUT:
            f = (scalar) frequency to evaluate at
            
        OUTPUT:
            WT[1,-1} = stress value at the earth's surface (r = rspan[1])
            
        UPDATES GLOBAL VARIABLES:
            omega = angular frequency for stress_disp_tor
            rvec = radii at which displacement and stress eigenfunctions were evaluated at
            WT = displacement and stress eigenfunctions"""
    
    #global omega, rvec, rspan, WT
    
    WT0 = np.array([1.0, 0.0]) # initial values of [displacement, stress]
    spshell_config.omega = 2*np.pi*f # angular frequency for stress_disp_tor
    
    # note: the dimension of rvec and WT is the number of points needed for
    # the numerical integration -- this will vary. You can adjust it via the 'max_step' parameter (try 1E2 or 1E3)
    rspan_t = tuple(spshell_config.rspan.tolist())
    sol = solve_ivp(stress_disp_tor,rspan_t,WT0,max_step=5E4)
    WT = sol.y
    rvec = sol.t
    spshell_config.WT = WT
    spshell_config.rvec = rvec
    
    return WT[1,-1]   # stress value at earth's surface (r = rspan[1])