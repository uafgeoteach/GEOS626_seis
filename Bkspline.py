import numpy as np

def Bkspline(clon, clat, q, lon_vec, lat_vec, ncol=1):
    """ INPUT:
            clon, clat, q   = these describe the local spherical spline basis function
            ncol            = how many columns of ff you want returned (derivatives)
                              (=1 default)
                              Can be 1, or 3,4,5
            lon_vec,lat_vec = datapoints at which you want the spherical spline evaluated

        OUTPUT:
            ff              = value of the spline function (and derivatives)
                              evaluated at the specified lon-lat points"""
        
    # convert to theta-phi
    deg = 180/np.pi
    ph     = clon/deg
    th     = (90-clat)/deg
    ph_vec = lon_vec/deg
    th_vec = (90-lat_vec)/deg
    
    # options and parameters -- q controls the scale (width) of the spline
    nf    = 2**q
    c72   = np.cos(72/deg)
    base  = np.arccos(c72 / (1 - c72))
    db    = base / nf
    zeps  = 1e-3*base     # determines whether a datapoint is ON a gridpoint
    
    # datapoint locations
    costh = np.cos(th_vec)
    sinth = np.sin(th_vec)
    ndata = len(th_vec)
    
    # r : delta/delta-bar in WD95
    delta = np.arccos( np.cos(th) * costh + np.sin(th) * sinth * np.cos(ph - ph_vec) )
    r   = delta/ db
    dif = r - 1
    
    # separate r into three parts: assign each element to one of four regions
    inds1 = np.flatnonzero(dif > 1)     # outside outer circle
    inds2 = np.flatnonzero((dif <= 1) & (dif >= 0))  # within outer ring
    inds3 = np.flatnonzero((dif > -1 + zeps) & (dif < 0))  # within inner circle
    inds4 = np.flatnonzero(dif <= -1 + zeps)   # ON the center point
    
    # check
    if len(inds1) + len(inds2) + len(inds3) + len(inds4) - len(dif) != 0:
        print(len(inds1))
        print(len(inds2))
        print(len(inds3))
        print(len(inds4))
        print(len(inds1)+ len(inds2) + len(inds3) + len(inds4))
        print(len(dif))
        raise Exception('Data points have not been partitioned correctly')
        
    if ncol == 1:
        ff = np.zeros((ndata,1))
        ff[inds2] = ((-0.25*dif[inds2]  + 0.75)*dif[inds2]  - 0.75)*dif[inds2] + 0.25
        ff[inds3] = (0.75*r[inds3] - 1.5) * r[inds3]**2  + 1
        ff[inds4] = 1
        
    else:
        cosdel = np.cos(th)*costh + np.sin(th) * sinth * np.cos(ph - ph_vec)
        sindel = np.sqrt(1 - cosdel*cosdel)
        cotdel = cosdel / sindel
        
        # ddelta/dphi and ddelta/dtheta (see MMA file wang_arc.nb)
        dadp = ( np.sin(th) * sinth * np.sin(ph_vec - ph) ) / sindel
        dadt = ( np.cos(th) * sinth - costh * np.sin(th) * np.cos(ph - ph_vec) ) / sindel
        
        # db : delta-bar in WD95
        # d_near varies for each gridpoint, due to irregularities in grids
        dq = 1 / db
        
        # columns of ff :
        # (1) f, function value
        # (2) df/dph
        # (3) df/dth
        # (4) surf_del2 -- depends only on delta
        # (5) |del f|   -- depends only on delta
    
        # datapoint is outside the outer circle
        ff = np.zeros((ndata,ncol))
        
        # datapoint is within the outer ring
        ff[inds2,0] = np.ravel(((-0.25*dif[inds2] + 0.75)*dif[inds2] - 0.75) * dif[inds2] + 0.25)
        ff[inds2,1] = np.ravel(dq * (-0.75 + 1.5*dif[inds2] - 0.75*dif[inds2]**2) * dadp[inds2])
        ff[inds2,2] = np.ravel(dq * (-0.75 + 1.5*dif[inds2] - 0.75*dif[inds2]**2) * dadt[inds2])
        
        if ncol >= 4:
            ff[inds2,3] = np.ravel(dq * (3 - 1.5*r[inds2] + cotdel[inds2] * (-0.75 + 1.5*dif[inds2] - 0.75*dif[inds2]**2)))
            ff[inds2,4] = np.ravel(0.75 * db**-3 * (2*db - delta[inds2])**2)
            
        # datapoint is within the inner circle
        ff[inds3,0] = np.ravel((0.75*r[inds3] - 1.5) * (r[inds3]**2) + 1)
        ff[inds3,1] = np.ravel(dq * (-3*r[inds3] + 2.25*r[inds3]**2) * dadp[inds3])
        ff[inds3,2] = np.ravel(dq * (-3*r[inds3] + 2.25*r[inds3]**2) * dadt[inds3])
        
        if ncol >= 4:
            ff[inds3,3] = np.ravel(dq * (-3 + 4.5*r[inds3] + cotdel[inds3] * (-3*r[inds3] + 2.25*r[inds3]**2)))
            ff[inds3,4] = np.ravel(0.75 * db**-3 * (4*db - 3*delta[inds3]) * delta[inds3])
            
        # datapoint is in the vicinity of the target spline centerpoint
        # FIX THIS : see Wang & Dahlen (1995)
        # here we simply assign it the closest value
        
        if len(inds4) > 0:
            if ncol > 3:
                igood = np.nonzero(dif > -1 + zeps)
                imin  = np.amin(r[igood])
                d2val = ff[imin,3]
                tvec = np.zeros((1,ncol))
                tvec[0] = 1
                tvec[-1] = d2val
                ff[inds4,0] = np.matlib.repmat(tvec,len(inds4),1)
                
            elif ncol == 3:
                ff[inds4,0:3] = np.array([1, 0, 0])
                
            elif ncol == 1:
                ff[inds4,0] = 1
                
    return ff