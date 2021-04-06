import numpy as np


def earthfun_love(r, rspan, crho, cmu, mrho, mmu, cthick):
    """Python adaptation of earthfun_love.m by Carl Tape.

    Amanda McPherson, Jan 2021
    Some edits by Liam Toney, April 2021

    earthfun_love() returns a rho and mu value for a specified radius
    r. Called by stress_disp_love().

    Args:
        r (int or float): Earth radius to compute rho and mu for
        rspan:
        crho (int or float):
        cmu (int or float):
        mrho (int or float):
        mmu (int or float):
        cthick (int or float):

    Returns:
        Tuple of (rho, mu)
    """

    # Enter code here!

    return rho, mu
