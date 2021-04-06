import numpy as np
from PREMfun import PREMfun


def earthfun(r, rspan, imod):
    """Python adaptation of earthfun.m by Carl Tape.

    Amanda McPherson, Dec 2020
    Some edits by Liam Toney, April 2021

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
        Prho = np.array([-9.9274*1E-04, 9.0147*1E+03])
        Pmu = np.array([-7.8035*1E+04, 5.6536*1E+11])
        rho = np.polyval(Prho,r)
        mu = np.polyval(Pmu,r)

    elif imod == 2:
        # cubic
        Prho = np.array([-2.84710*1E-16, 3.84976*1E-9, -1.76479*1E-2, 3.24479*1E4])
        Pmu  = np.array([-8.11871*1E-9, 9.56717*1E-2, -4.250608*1E5, 9.578569*1E11])
        rho = np.polyval(Prho,r)
        mu = np.polyval(Pmu,r)

    elif imod == 3:
        #if r is in the water, return the uppermost crustal value
        r0 = np.array(r)
        r0[r0 >= rspan[1] - 3000] = rspan[1] - 3001
        #rho = VelocityModel.evaluate_below(r0,'r')
        #Svel = VelocityModel.evaluate_below(r0,'s')
        #mu = rho*Svel*Svel
        # Only need rho and mu from this
        iz,a,b,r,Qm,Qk,acorr,bcorr,k,m,Qa,Qb,rbound = PREMfun(rspan[1]-r0)
        rho = r
        mu = m

    else:
        raise ValueError('imod must be one of [1, 2, 3]')

    return rho, mu
