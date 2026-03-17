import numpy as np

def earthfun(r, rspan, imodel):
    """
    earthfun() returns a rho and mu value for a specified radius r (and rspan
    and imodel). Called by stress_disp_tor().

    Args:
        r (int or float): Earth radius to compute rho and mu for
        rspan:
        imodel (int):

    Returns:
        Tuple of (rho, mu)
    """

    if imodel == 1:
        # linear
        # enter your code here
        cmbr = rspan[0]    # b
        earthr = rspan[1]  # a
        
        raise Exception('earthfun.py imodel=1 not yet implemented (comment this out when you have implemented it)')

    elif imodel == 2:
        # cubic
        # enter your code here
        
        raise Exception('earthfun.py imodel=2 not yet implemented (comment this out when you have implemented it)')

    elif imodel != 1 or imodel != 2:
        raise Exception('invalid imodel (=1,2)')

    return rho, mu
