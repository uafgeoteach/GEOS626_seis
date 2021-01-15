import numpy as np
from scipy.integrate import solve_ivp
from stress_disp_love import stress_disp_love
import spshell_love_config

def surf_stress_love(k0):
    """ Modified from surf_stress_love.m by Carl Tape
        Adapted by Amanda McPherson (Jan 2021)
        SURF_STRESS_LOVE calculate stress at the surface of a toroidal earth
        Applied Seismology, GEOS 626, University of Alaska Fairbanks

        modified from surf_stress.py in the modes homework
        
        INPUT:
            k0 - wavenumber, scalar
            
        OUTPUT:
            WT[1,-1] - stress value at earth's surface (r = rspan(2)), scalar
            
        UPDATES:
            spshell_love_config.k
            spshell_love_config.rvec
            spshell_love_config.WT"""
    
    spshell_love_config.k = k0
    k = spshell_love_config.k
    omega = spshell_love_config.omega
    mbeta = spshell_love_config.mbeta
    mmu = spshell_love_config.mmu
    
    # calculate initial conditions at r=0 within the mantle halfspace
    mk = omega/mbeta
    nub = np.sqrt(k**2 - mk**2)
    
    if np.iscomplex(nub):
        print('setting nub=0 (k=%.3e mk=%.3e)'% (k,mk))
        nub = 0
        
    Tbot = mmu*nub
    WT0 = np.array([1.0, Tbot])   # the initial values of [displacement stress]
    
    # note: the dimension of rvec and WT is the number of points needed for
    # the numerical integration -- this will vary. You can adjust it via the 'max_step' parameter
    rspan_t = tuple(spshell_love_config.rspan.tolist())
    sol = solve_ivp(stress_disp_love,rspan_t,WT0,max_step=5.3*1E3)
    WT = sol.y
    rvec = sol.t
    spshell_love_config.WT = WT
    spshell_love_config.rvec = rvec
    
    return WT[1,-1]   # stress value at earth's surface (r = rspan[1])
        