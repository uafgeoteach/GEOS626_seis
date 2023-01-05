import numpy as np

from earthfun import earthfun

def stress_disp_tor(r, WT_0, l, omega, imod, rho, mu, rspan):
    """Python adaptation of stress_disp_tor.m by Carl Tape.

    Coding by Amanda McPherson, Dec 2020
    Some edits by Liam Toney, April 2021

    Finds the derivatives of W(r) and T(r)

    Args:
        r (int or float): Radius at which to evaluate the derivatives
        WT_0 (numpy.ndarray): WT_0[0] is W(r) and WT_0[1] is T(r)
        l (int):
        omega (int or float): Angular frequency
        imod (int):
        rho:
        mu:
        rspan:

    Returns:
        numpy.ndarray: dWT. dWT[0] is the derivative of W(r), and dWT[1] is the derivative of T(r)
    """

    dWT = np.empty(2)
    # structural values at radius r: density and rigidity
    # note: if imod == 0, then we use the rho and mu from spshell.ipynb
    if imod != 0:
        rho, mu = earthfun(r, rspan, imod)

    # displacement (first row of equation 1)
    dWT[0] = WT_0[0] / r + WT_0[1] / mu

    # stress (second row of equation 2)
    dWT[1] = ((l-1)*(l+2)*mu/(r*r) - rho*omega*omega)*WT_0[0] - 3*WT_0[1]/r

    return dWT
