
import numpy as np
import math as m

# python version of Matlab's sph2cart and cart2sph function
# transforms elements of the spherical coordinate arrays azimuth, elevation, and r to Cartesian, or xyz, coordinates.

def sph2cart(azimuth,elevation,r):
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z
