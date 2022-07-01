import numpy as np
import matplotlib.pyplot as plt

def plot_histo(hdat,edges,itype=2,make_plot=True):
    #PLOT_HISTO plot a histogram with cyan bars and black boundaries
    #
    # INPUT
    #   hdat        input data to bin
    #   edges       vector defining the edges of the bins (for hdat)
    #   itype       optional: type of histogram (=1,2,3) [default = 2]
    #   make_plot   optional: plot histogram [default = true]
    #hdat = hdat.flatten();
    #barcolor = [1, 1, 1]*0.8;
    
    # bin width (only relevant if bins are the same width)
    dbin = edges[1] - edges[0]
    hedges=np.append(edges,edges[-1]+5)
    Ntotal = len(hdat);
    N,b = np.histogram(hdat,hedges);
    if itype ==1:
        Nplot = N; xlab = 'Count'
    if itype ==2: 
        Nplot = np.divide(N,Ntotal); xlab = 'Fraction'
    if itype ==3: 
        Nplot = np.divide(np.divide(N,Ntotal),dbin); xlab = 'PDF'
        #if len(unique(edges)) > 1:
        if np.std(np.diff(edges))/np.mean(np.diff(edges)) > 1e-4:       # ad hoc criterion
            print(np.unique(np.diff(edges)))
            print('PDF is not implemented to allow bins with varying widths')
            
    elif itype!=1 and itype!=2 and itype!=3: 
        print('itype = %i -- it must be 1,2, or 3'%(itype)) 

    if make_plot==True:
        plt.bar(edges,Nplot, width=0.8*dbin);
        plt.xlim([min(edges), max(edges)]);
        plt.ylabel('%s (N=%i)'% (xlab,Ntotal))
        
        if len(hdat) != np.sum(N):
            print('(plot_histo.m): You may want to extend the histogram edges -->');
            print(' there are %i/%i input that are outside the specified range'%
                (len(hdat)-np.sum(N),len(hdat)))
            #disp(sprintf(' the number of input (%i) does not equal the sum of bin counts (%i).',length(hdat),sum(N)));
    
    plt.tight_layout()