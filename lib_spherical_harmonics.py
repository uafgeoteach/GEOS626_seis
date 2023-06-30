import numpy as np
import matplotlib.pyplot as plt

from scipy.special import lpmv
from scipy.interpolate import griddata

###############################################################################################################

def griddataX(phi, lat, alphcent, npts):
    """
    This function takes lat-lon points (radians) and a function
    defined at each point, and it returns a set of points
    XYZ that can be plotted as a planar surface plot.
    Note that the lat-range is 180 and the lon-range is 360.

    phi       = longitude (radians)
    lat       = latitude (radians)
    alphcent  = function values defined at every grid point
    npts      = number of points to use in horizontal interpolation

    griddataX.m by Carl Tape
    Python adaptation by Amanda McPherson
    """

    deg = 180 / np.pi

    # interpolation commands for surface/contour plots
    x = phi * deg
    y = lat * deg
    XY = np.append(x, y, axis=1)
    z = alphcent  # z = alphcent.flatten()

    xlin = np.linspace(min(x), max(x), npts)
    xlin = np.reshape(xlin, (len(xlin), 1))
    ylin = np.linspace(min(y), max(y), npts)
    ylin = np.reshape(ylin, (len(ylin), 1))

    X, Y = np.meshgrid(xlin, ylin)
    Z = griddata(XY, z, (X, Y), method='cubic', fill_value=0)

    return X, Y, Z

###############################################################################################################

def xlm(th, LL):
    """
    XLM generalized Legendre functions (xlm), Dahlen and Tromp p. 848
    Optional plots of L+1 xlm (0 <= M <= L), which requires you
    to use ifigs = 1 and execute as a script.

    Note: to look at legendre functions, use legplots3.m

    th = vector of polar angles (0 to pi)
    LL = degree for the spherical harmonic

    see Tape (2003) appendix in thesis for ylm notation and conventions

    xlm.m by Carl Tape, 2006-01-22

    Python adaptation by Amanda McPherson, 2021
    """

    num = len(th)
    X = np.cos(th)

    m = np.arange(LL + 1)
    # THIS IS THE PYTHON CALCULATION OF THE ASSOCIATED LEGENDRE FUNCTIONS.
    # It calculates the functions from M = 0 to M = L.
    # We have to modify it to be consistent with the useage in JHW's programs:
    #       1. Multiply everything by factor A.
    #       2. Change the sign of every other order.
    #       3. Divide all non-M=0 values by sqrt(2)

    A = np.sqrt((2 * LL + 1) / (4 * np.pi))  # normalization (amplitude correction)
    legmat = lpmv(m, LL, X)
    legmat = A * legmat

    for MM in range(LL + 1):
        legmat[:, MM] = (-1) ** MM * legmat[:, MM]  # changes the sign of every other column of values

        if MM != 0:
            norm = ((-1) ** MM) * np.sqrt(
                (2 * np.math.factorial(LL - MM)) / np.math.factorial(LL + MM))  # divides by sqrt(2) if M is not 0
            legmat[:, MM] = norm * legmat[:, MM] / np.sqrt(2)  # Applies Schmidt normalization

    return legmat

###############################################################################################################

def ylm(th, phi, L, M):
    """
    YLM surface spherical harmonic functions Y_lm

    Convention similar to Dahlen and Tromp (1999), p. 851.
    th    = vector of polar angles [0:pi]
    phi   = vector of azimuthal angles [-pi:pi]
    L     = degree
    M     = order [-L <= M <= L]

    calls xlm

    ylm.m by Carl Tape, 2003-02-24
    Python adaptation by Amanda McPherson, 2021
    """

    num = len(th)

    if M == 0 and L == 0:
        raise Exception('L=0, M=0 is the breathing mode, which is spherically symmetric')

    pts = xlm(th, L)  # key command
    pts2 = pts[:, M]  # python has M in the columns already

    # automatically does all the orders (0 <= M <= L), returns a flattened ndarray of length num
    if M == 0:
        val = pts2

    elif M > 0 and M <= L:
        val = pts2 * np.sin(M * phi).flatten()

    else:
        raise Exception('invalid M or L: requirement that 0 <= M <= L')

    return val

###############################################################################################################

def ylmplots_fun():
    """
    Plots Ylm real spherical harmonic functions, either as a 3D figure or a 2D contour.
    ylmplots.m by Carl Tape, 2004-04-13
    Python adaptation by Amanda McPherson
    calls griddataX
    """

    stit = 'L=%i, M=%i : (%i, %i)' % (L, M, 2 * M, L - M + 1)
    lat = np.pi / 2 - th

    X, Y, Z = griddataX(ph, lat, alphcent, npts)
    zmax = np.amax(abs(Z))

    plt.contour(X, Y, Z, levels=numc, vmin=-zmax, vmax=zmax)
    plt.xticks([-180, 0, 180])
    plt.yticks(thtick)
    plt.title(stit)
    plt.xlabel('Longitude (φ)')
    plt.ylabel('Latitude (90 - θ)')
    plt.axis([-180, 180, -90, 90])
    plt.colorbar()

###############################################################################################################