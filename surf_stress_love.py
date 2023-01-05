import numpy as np

from scipy.integrate import solve_ivp

from stress_disp_love import stress_disp_love

def surf_stress_love(k0, omega, rspan, crho, cmu, mrho, mmu, cthick, max_step=5e3, return_wt_rvec=False):
    """
    surf_stress_love() calculates the stress at the surface of a layer-over-halfspace

    Args:
        k0 (int or float): Wavenumber
        omega:
        mmu:
        rspan:
        crho:
        cmu:
        mrho:
        cthick:
        max_step (int or float): Step size for ODE solver
        return_wt_rvec (bool): If True, return WT and rvec in addition to WT[1,-1]

    Returns:
        WT[1,-1] = stress value at the Earth's surface (r = rspan[1])
        [WT = displacement and stress eigenfunctions]
        [rvec = radii at which displacement and stress eigenfunctions were evaluated at]
    
    Applied Seismology (GEOS 626) University of Alaska Fairbanks
    contributors: Charles Ammon, Carl Tape, Amanda McPherson*, Liam Toney
    """

    k = k0

    # calculate initial conditions at r=0 within the mantle halfspace
    mbeta = np.sqrt(mmu/mrho)
    mk    = omega/mbeta
    nub   = np.sqrt(k**2 - mk**2)

    if np.iscomplex(nub):
        print('setting nub=0 (k=%.3e mk=%.3e)'% (k,mk))
        nub = 0

    Tbot = mmu*nub
    WT0  = np.array([1.0, Tbot])   # the initial values of [displacement stress]

    # note: the dimension of rvec and WT is the number of points needed for
    # the numerical integration -- this will vary. You can adjust it via the 'max_step' parameter
    rspan_t = tuple(rspan.tolist())
    sol = solve_ivp(stress_disp_love, rspan_t, WT0, max_step=max_step, args=(k, omega, rspan, crho, cmu, mrho, mmu, cthick))
    WT = sol.y
    rvec = sol.t

    if return_wt_rvec:
        return WT[1,-1], WT, rvec
    else:
        return WT[1,-1]   # stress value at Earth's surface (r = rspan[1])
