import numpy as np

from scipy.integrate import solve_ivp

from stress_disp_tor import stress_disp_tor

def surf_stress(f, l, imod, rho, mu, rspan, max_step=5e4, return_wt_rvec=False):
    """
    Integrates using the RK45 method, and calls stress_disp_tor() to calculate the
    derivatives of displacement and stress

    Args:
        f (int or float): Frequency to evaluate at
        l:
        imod:
        rho:
        mu:
        rspan:
        max_step (int or float): Step size for ODE solver
        return_wt_rvec (bool): If True, return WT and rvec in addition to WT[1,-1]

    Returns:
        WT[1,-1] = stress value at the Earth's surface (r = rspan[1])
        [WT = displacement and stress eigenfunctions]
        [rvec = radii at which displacement and stress eigenfunctions were evaluated at]
    """

    WT0 = np.array([1.0, 0.0]) # initial values of [displacement, stress]
    omega = 2*np.pi*f # angular frequency for stress_disp_tor()

    # note: the dimension of rvec and WT is the number of points needed for
    # the numerical integration -- this will vary. You can adjust it via the 'max_step' parameter (try 1e2 or 1e3)
    rspan_t = tuple(rspan.tolist())
    sol = solve_ivp(stress_disp_tor,rspan_t,WT0,max_step=max_step, args=(l, omega, imod, rho, mu, rspan))
    WT = sol.y
    rvec = sol.t

    if return_wt_rvec:
        return WT[1,-1], WT, rvec
    else:
        return WT[1,-1]   # stress value at Earth's surface (r = rspan[1])
