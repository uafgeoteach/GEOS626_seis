import numpy as np

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