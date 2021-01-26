import numpy as np

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