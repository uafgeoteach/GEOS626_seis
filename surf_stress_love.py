import numpy as np
from scipy.integrate import solve_ivp
from stress_disp_love import stress_disp_love


def surf_stress_love(k0, omega, mbeta, mmu, rspan, crho, cmu, mrho, cthick, max_step=5.3e3, return_wt_rvec=False):
    """Modified from surf_stress_love.m by Carl Tape.

    Adapted by Amanda McPherson (Jan 2021)
    Some edits by Liam Toney, April 2021

    surf_stress_love() calculates the stress at the surface of a
    toroidal Earth.

    Applied Seismology (GEOS 626) University of Alaska Fairbanks

    Modified from surf_stress.py in the modesA homework

    Args:
        k0 (int or float): Wavenumber
        omega:
        mbeta:
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
    """

    k = k0

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
    rspan_t = tuple(rspan.tolist())
    sol = solve_ivp(stress_disp_love,rspan_t,WT0,max_step=max_step, args=(k, omega, rspan, crho, cmu, mrho, mmu, cthick))
    WT = sol.y
    rvec = sol.t

    if return_wt_rvec:
        return WT[1,-1], WT, rvec
    else:
        return WT[1,-1]   # stress value at Earth's surface (r = rspan[1])
