#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# Library of functions for the GEOS 626 course at UAF
# Contributers: Carl Tape, Nealey Sims*, Amanda McPherson*, Aakash Gupta*

import datetime as dt
import matplotlib.pyplot as plt
import math as m
import numpy as np
import scipy.io

from itertools import product, combinations

############################################################

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
    
    # allow for scalar input values
    lon_vec = np.atleast_1d(lon_vec)
    lat_vec = np.atleast_1d(lat_vec)
    
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
        if len(inds2) != 0:
            ff[inds2] = ((-0.25*dif[inds2]  + 0.75)*dif[inds2]  - 0.75)*dif[inds2] + 0.25
        
        if len(inds3) != 0:
            ff[inds3] = (0.75*r[inds3] - 1.5) * r[inds3]**2  + 1
            
        if len(inds4) != 0:
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

############################################################

def globefun3(R,lat,lon,bool_point,lc, fig,ax):
    # Python adaptation of globefun3.m
    # Python coding done by Nealey Sims
    #
    # GLOBEFUN3 plots a single point on the sphere with a reference globe
    #
    # INPUT
    #   R           radius of sphere
    #   lat         latitude (deg)  (theta = 90 - lat)
    #   lon         longitude (deg) (lon = phi)
    #   bool_point  plot the specified point
    #   lc          line style/color
    #
    #
    
    deg = 180/np.pi
    lw = 1.0
    th = (90-lat)/deg
    ph = lon/deg
    
    # (1) xy-plane circle (equator)
    # (2) yz-plane
    # (3) xz-plane (prime meridian)
    gpts = 100
    tt = np.linspace(0, 2*np.pi, gpts)
    
    
    ax.set_box_aspect((4,4,3.4))  ## Equal axis
    ax.plot(R*np.cos(tt), R*np.sin(tt), np.zeros(gpts), lc, linewidth=lw)
    ax.plot(np.zeros(gpts), R*np.cos(tt), R*np.sin(tt), lc, linewidth=lw)
    ax.plot(R*np.cos(tt), np.zeros(gpts), R*np.sin(tt), lc, linewidth=lw)
    
    fac = 1.15; lw2 = 0.5
    RR = np.array([-R,R])
    ax.plot((fac*RR),[0, 0],[0, 0],'k', linewidth=lw2)
    ax.plot([0, 0],(fac*RR),[0, 0],'k', linewidth=lw2)
    ax.plot([0, 0],[0, 0],(fac*RR),'k', linewidth=lw2)
    
    # draw cube
    r = (fac*RR)
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), color="k", linewidth=lw2)
    
    
    if bool_point == 1:
        # plot the specified point P(r,th,phi)
        msize = 20; lw3 = lw+0.5
        # [X,Y,Z] = sph2cart(az,el,r)
        [X,Y,Z] = sph2cart(ph,np.pi/2-th,R)  #Note: need sph2cart function in folder
        
        [X2,Y2,Z2] = sph2cart(ph,0,R)
  
        ax.plot([X],[Y],[Z], color=lc, marker='o')
        ax.plot([0, X],[0, Y],[0, Z],lc, linewidth= lw3)
        
        leps = 0.1;
        if ( lat >= -90+leps) and (lat <= 90-leps):
            xzplane = np.array([R*np.cos(tt) , np.zeros(gpts) , R*np.sin(tt)])
            rotmat = np.array([[np.cos(ph), -np.sin(ph), 0] , [np.sin(ph), np.cos(ph), 0], [0, 0, 1]])
            rotxz = np.dot(rotmat,xzplane)
            ax.plot(rotxz[0,:], rotxz[1,:], rotxz[2,:],lc,linewidth=lw);
            aa = lat/deg;
            ax.plot(R*np.cos(aa)*np.cos(tt), R*np.cos(aa)*np.sin(tt), R*np.sin(aa)*np.ones(gpts),lc,linewidth =lw);
            
            # plots some extra lines for reference
            ax.plot([0, X2],[0, Y2],[0, 0],lc, linewidth= lw2);
            ax.plot([X, X],[Y, Y],[0, Z],lc, linewidth= lw2);
            ax.plot([0, X],[0, Y],[Z, Z],lc, linewidth= lw2);
        #plt.title('\u03C6 = %.2f, \u03B8 = %.2f' % (ph*deg, th*deg))
        ax.view_init(10, 10)
    ax.grid(False)
    ax.xaxis.set_pane_color((0.0, .0, 0.0, 0.0))
    ax.yaxis.set_pane_color((0.0, .0, 0.0, 0.0))
    ax.zaxis.set_pane_color((0.0, .0, 0.0, 0.0))
    #plt.show()

############################################################

def matlab2datetime(matlab_datenum):
        # equivalent of datestr when converting from serial number to date
       
        day = dt.datetime.fromordinal(int(matlab_datenum))
        dayfrac = dt.timedelta(days=float(matlab_datenum)%1) - dt.timedelta(days = 366)
        return day + dayfrac
    
############################################################

# To use function: cid=fig.canvas.mpl_connect('button_press_event', markp)
# 
# This marks the 1/x value on a plot where you click the mouse.
# python version based on the markt.m writtin by Carl Tape
#
# INPUT ARGUMENTS:
# 
# left mouse button click: plot point
#
# EXAMPLE:
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlim([0, 10])
# ax.set_ylim([0, 10]) 
# markt(fig)
#
# pt coordinates will print to screen
  
def markp(event):
    print('x=%f, y=%f' % (event.xdata, event.ydata))
    prd=round(1/event.xdata, 3)
    axe=event.inaxes
    axe.text(event.xdata, event.ydata, s = str(prd))
    #plt.plot(event.xdata, event.ydata, 'ro')
    plt.draw()    
    
def markp_min(event):
    print('x=%f, y=%f' % (event.xdata, event.ydata))
    prd=round(1/event.xdata/60, 1)
    axe=event.inaxes
    axe.text(event.xdata, event.ydata, s = str(prd))
    #plt.plot(event.xdata, event.ydata, 'ro')
    plt.draw()    
    
############################################################

def seis2GR(mag, dmag, idisplay=1, ifigure=0):
    """ Python version of seis2GR.m by Carl Tape.
        Converts seismicity catalog to Gutenberg-Richter frequency-magnitude distribution, in both a cumulative
        and not cumulative version.
        
        INPUT:
            mag - ndarray of magnitudes for a set of events
            dmag - magnitude increment for bins in histogram
        OPTIONAL:
            idisplay - display numbers in each bin (default=1)
            ifigure - plot different versions of histogram (default=0)(See also GR2plot.m(THIS SHOULD BE IN PYTHON EVENTUALLY))
        
        OUTPUT:
            N - cumulative number of events (from largest to smallest magnitude)
            Ninc - number of events per magnitude bin
            Medges - edges of the magnitude bins
        
        Dependencies:
            numpy
            scipy.io
            matplotlib.pyplot"""
    
    n = len(mag)
    minm = min(mag)
    maxm = max(mag)
    
    print('seis2GR: %i events, min M = %.3f, max M = %.3f'% (n,minm,maxm))
    
    emin = np.floor(minm/dmag)*dmag
    emax = (np.ceil(maxm/dmag)+1)*dmag # +1 in case maxm is right on an edge
    Medges = np.arange(emin,emax,dmag).tolist()
    Medges = np.round(Medges, decimals=2)
    
    # see numpy documentation for details: the bins are closed on the left and open on the right
    Ninc, bin_edges = np.histogram(mag,Medges)
    N = np.flipud(np.flipud(Ninc).cumsum())
    nbin = len(Ninc)
    
    if idisplay == 1:
        for ii in range(nbin): #range(x) goes from 0 to x-1
            print('bin ',ii,': Mw = [',Medges[ii],' ',Medges[ii+1],'] Ninc = ',Ninc[ii],'N = ',N[ii])
    
    if ifigure == 1:
        plt.figure(figsize=(9,8))
        for kk in [1, 2, 3, 4]:
            if kk == 1:
                D = Ninc
                ylab = 'Number (N = '+str(n)+')'
            if kk == 2:
                D = np.log10(Ninc)
                D[np.isinf(D)] = 0 #plt.hist does not work if any value is +/-inf
                ylab = 'Log10[Number] (N = '+str(n)+')'
            if kk == 3:
                D = N
                ylab = 'Cumulative Number (N = '+str(n)+')'
            if kk == 4:
                D = np.log10(N)
                D[np.isinf(D)] = 0
                ylab = 'Log10[Cumulative Number] (N = '+str(n)+')'
            plt.subplot(2,2,kk)
            plt.hist(Medges[:-1], bins=Medges, weights=D)
            plt.xlabel('Magnitude')
            plt.ylabel(ylab)
        plt.show()
        
    return N, Ninc, Medges

############################################################

def smooth(a,WSZ):
    """ Python adaptation of MATLAB smooth.m function
        Author: Divakar on Stack Overflow
        
        INPUT:
            a: Numpy 1-D array containing the data to be smoothed
            WSZ: smoothing window size needs. Must be an odd number, as in the original implementation
            
        OUTPUT:
            Numpy 1-D array of smoothed data points"""
    
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    
    return np.concatenate((  start , out0, stop  ))

############################################################

# python version of Matlab's sph2cart and cart2sph function
# transforms elements of the spherical coordinate arrays azimuth, elevation, and r to Cartesian, or xyz, coordinates.

def sph2cart(azimuth,elevation,r):
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

############################################################

def wf_fft(wf,fNyq):
    """ Python adaptation of wf_fft.m by Michael West
        Necessary for GEOS626 work
        
        INPUT:
            wf - Numpy array of the data points in your trace
            fNyq - the Nyquist frequency
        
        OUTPUT:
            fft_amp - Numpy array of spectral amplitudes
            fft_phase - Numpy array of phases
            fft_freq - Numpy array of frequencies"""
    
    NFFT = int(2**(np.ceil(np.log(len(wf))/np.log(2))))  # Next highest power of 2
    FFTX = np.fft.fft(wf,n=NFFT)                       # Take fft, padding with zeros.
    NumUniquePts = int(np.ceil((NFFT+1)/2))
    FFTX = FFTX[0:NumUniquePts]              # throw out neg frequencies
    MX = abs(FFTX)                           # Take magnitude of X
    MX = MX*2                                # Multiply by 2 
    fft_amp = MX/len(wf) 
    
    fft_phase = np.angle(FFTX)               # Take magnitude of X
    
    f = (np.arange(NumUniquePts))*2/NFFT            
    fft_freq = f*fNyq
    
    return fft_amp, fft_phase, fft_freq

############################################################
