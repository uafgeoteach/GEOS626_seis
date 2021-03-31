import numpy as np

def gridvec(xmin,xmax,numx,ymin,ymax):
    """  This function inputs specifications for creating a grid 
         of uniformly spaced points, reshaped into column vectors
         of the x- and y-coordinates.  Note that dx = dy."""
    
    xvec0 = np.linspace(xmin,xmax,numx)
    dx = xvec0[1] - xvec0[0]
    yvec0 = np.arange(ymin, ymax+dx, dx)
    
    X, Y = np.meshgrid(xvec0,yvec0)
    a,b = X.shape
    xvec = np.reshape(X,(a*b,1))
    yvec = np.reshape(Y,(a*b,1))
    
    numy = len(yvec0)
    
    return xvec, yvec, numy, X, Y

#def gridvec(xmin,xmax,numx,ymin,ymax,returnXY=False):
#    """  This function inputs specifications for creating a grid 
#         of uniformly spaced points, reshaped into column vectors
#         of the x- and y-coordinates.  Note that dx = dy."""
#    
#    xvec0 = np.linspace(xmin,xmax,numx)
#    dx = xvec0[1] - xvec0[0]
#    yvec0 = np.arange(ymin, ymax+dx, dx)
#    
#    X, Y = np.meshgrid(xvec0,yvec0)
#    a,b = X.shape
#    xvec = np.reshape(X,(a*b,1))
#    yvec = np.reshape(Y,(a*b,1))
#    
#    numy = len(yvec0)
#    
#    if returnXY:
#        return xvec, yvec, numy, X, Y
#    else:
#        return xvec, yvec, numy