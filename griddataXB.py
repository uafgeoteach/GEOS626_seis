import numpy as np
from scipy.interpolate import griddata

def griddataXB(xvec,yvec,zvec,npts,stype='linear'):
    """ Converts irregular points to regular points, for mesh plotting.

            INPUT:
                xvec,yvec   = map coordinates of irregular points
                zvec        = function values at (xi,yi)
                npts        = number of points in horizontal interpolation
                stype       = type of interpolation

                    'linear'    - Triangle-based linear interpolation (default).
                    'cubic'     - Triangle-based cubic interpolation.
                    'nearest'   - Nearest neighbor interpolation.

            OUTPUT:
                X,Y         = interpolated mesh
                Z           = interpolated function """
    
    # Reshape map coords for griddata
    xy_vec = np.append(xvec,yvec,axis=1)
    
    # construct mesh with UNIFORM spacing in x and y directions
    xlin  = np.linspace(min(xvec), max(xvec), npts)
    dx    = xlin[1] - xlin[0]
    ylin  = np.arange(min(yvec), max(yvec)+dx, dx)
    X, Y = np.meshgrid(xlin,ylin)
    
    # determine interpolated function using xvec,yvec input
    Z = griddata(xy_vec,zvec,(X,Y),method=stype,fill_value=0)
    
    return X, Y, Z