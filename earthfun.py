import numpy as np

def earthfun(r, rspan, imod):
    """
    earthfun() returns a rho and mu value for a specified radius r (and rspan
    and imod). Called by stress_disp_tor().

    Args:
        r (int or float): Earth radius to compute rho and mu for
        rspan:
        imod (int):

    Returns:
        Tuple of (rho, mu)
    """

    if imod == 1:
        # linear
        # enter your code here
        cmbr = rspan[0]    # b
        earthr = rspan[1]  # a
        
        raise Exception('earthfun.py imod=1 not yet implemented (comment this out when you have implemented it)')

    elif imod == 2:
        # cubic
        # enter your code here
        
        raise Exception('earthfun.py imod=2 not yet implemented (comment this out when you have implemented it)')

    elif imod != 1 or imod != 2:
        raise Exception('invalid imod (=1,2)')

    return rho, mu
