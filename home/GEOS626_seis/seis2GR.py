import numpy as np
import matplotlib.pyplot as plt
import scipy.io

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
    
    Ninc, bin_edges = np.histogram(mag,Medges)
    N = np.flipud(np.flipud(Ninc).cumsum())
    nbin = len(Ninc)
    
    if idisplay == 1:
        for ii in range(nbin): #range(x) goes from 0 to x-1
            print('bin ',ii,': Mw = [',Medges[ii],' ',Medges[ii+1],'] Ninc = ',Ninc[ii],'N = ',N[ii])
    
    if ifigure == 1:
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
                ylab = 'Cumulative[Number] (N = '+str(n)+')'
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
