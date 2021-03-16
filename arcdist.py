import numpy as np

def arcdist(slon,slat,rlon,rlat):
    """ This function determines the arc distance from an earthquake source (s) to a seismic receiver (r),
         IGNORING THE OBLATENESS OF THE EARTH (see arcd.m)
         Based on Matlab function of the same name. 
         
         INPUT (all in DEGREES):
            slat = latitude of source
            slon = longitude of source
            rlat = latitude of receiver (ndarray)
            rlon = longitude of receiver (ndarray)

         OUTPUT (in degrees):
            dvec = arc-distances from (slat, slon) to each (rlat, rlon) """
    
    deg = 180/np.pi
    nrec = len(rlat)
    nsrc = len(slat)
    
    if nsrc == 1:
        slat = slat * np.ones(nrec)
        slon = slon * np.ones(nrec)
    else:
        raise Exception('Assume you only have one source, i.e. slat and slon are scalars')
        
    las = slat / deg
    lar = rlat / deg
    lod = (slon-rlon)/deg 
    
    dvec = deg * np.arccos( np.sin(lar)*np.sin(las) + np.cos(lar)*np.cos(las)*np.cos(lod) )
    
    return dvec