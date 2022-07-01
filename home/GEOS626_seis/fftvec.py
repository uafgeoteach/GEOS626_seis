import numpy as np

def fftvec(t):
    # Python version of fftvec.m made by Carl tape
    # FFTVEC provides frequency vector for Matlab's fft convention
    # 
    # This is the bare-bones FFT example for teaching purposes in the sense
    # that there is no padding with zeros, no odd-number time series, no
    # all-positive frequencies, and no use of fftshift.
    #
    # EXAMPLE:
    #   t = np.arange(2,9.0,0.4); n=len(t); f = fftvec(t); ix = 4; 
    #   dt = t[1]-t[0]; df = 1/(2*dt)/(n/2) 
    #   fig=plt.figure();plt.plot(f,'-'); 
    #   fig2=plt.figure(); plt.plot(np.fft.fftshift(f),'-');
    #
    # NOTE: I chose to set fNyq as the negative frequency such that
    #       fftshift(fftvec(t)) will give a vector increasing uniformly from
    #       -fNyq to just less than fNyq, with f=0 at index n/2+1
    #
    # NOTE: If you want to consider all-positive frequencies, then use abs(f)
    # 
    n = len(t)
    if n % 2 == 1:
       print('error(time series must have an even number of points)')

    dt = t[2] - t[1];       # sample rate
    fNyq = 1/(2*dt);        # Nyquist frequency
    
    # first half of frequency vector (note: first entry is f=0)
    f1 = np.transpose(np.linspace(0, float(fNyq), int((n/2)+1)));
    # full frequency vector
    #f = [f1 ; -f1(end-1:-1:2)];        % fNyq > 0
    
    f = np.concatenate([f1[0:int(n/2)] , -f1[:0:-1]])    # fNyq < 0
    return f
    
    # alternatively (for fNyq < 0)
    #df = 1/(n*dt);          % =2*fNyq/n
    #f1 = linspace(-fNyq,fNyq-df,n)'
    #f = [f1(n/2+1:n) ; f1(1:n/2)];

    #==========================================================================
    

