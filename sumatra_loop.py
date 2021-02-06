import matplotlib.pylab as plt
import numpy as np
import statistics as stats
from obspy.signal.invsim import cosine_taper
from matplotlib.mlab import detrend
from obspy import *
import os
from os import path
from wf_fft import wf_fft

def sumatra_loop(w0):
    # sumatra_loop.py
    #
    # Extract the time series for the 2004 Sumatra earthquake that are not
    # corrupted by obvious problems over the time window.
    #
    # This is called by sumatra_modes. The reason for separating it is
    # that it avoids having to reload the waveform data each time you run the
    # loop.
    #
    # To make composite PDF:
    # for file in `ls sumatra_modes*.ps` ; do ps2pdf $file ; done ; pdcat -r sumatra_modes*pdf all_sumatra_modes.pdf
    #
    
    savefile=False           ### change to False if sumatra_modes.txt file already saved
    compute_fft=True         ### create folder of computed fft
    spdy = 86400
    w = w0.copy()
    nw=len(w)
    print(nw)
    stas=[]
    chans=[]
    nets=[]
    locs=[]
    for tr in w:
        stas.append(tr.stats.station)
        chans.append(tr.stats.channel)
        locs.append(tr.stats.location)
        nets.append(tr.stats.network)
        tr.trim(tr.stats.starttime+(9.5*spdy),tr.stats.starttime+(19*spdy), pad=True, fill_value=0)
    # write to command window
    for ii in range (len(w)):
        print('%3i %7s %7s %4s %4s' % (ii,stas[ii],chans[ii],locs[ii],nets[ii]))
    
    
    # plot time series
    ifigure = 0
    
        # fill gaps in the time series if the total length of gaps is leff than
    # this fraction of the full time series.
    FTHRESH = 0.03
    
    imin = 0; imax = nw
    ii=imin
    scut=[]
    
    while ii < imax:
        figname='fig'+str(ii)
        stag = str(stas[ii])+ '_' +str(chans[ii])+ '_' +str(nets[ii])
        stdur = ('duration = %.2f days'+'duration') ### get w(ii)
        print('---------------------------------------')
        print('%i/%i %s' % (ii,nw,stag))
        stit = [str(stag) +', ' +str(stdur)]
    
        nd = len(w[ii].data)
    
        if ifigure==1: 
            fig=plt.figure()
            plt.plot(w[ii].data) 
            plt.title(stit)
    
        # specific commands for certain records
        
        if str(stas[ii])=='QSPA' and str(chans[ii])=='LHZ' and str(locs[ii]) == '20':
            print('special cutting of the end of the record')
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+4.48*1e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='ABKT':
            print('special cutting of the end of the record')
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+6.741e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
       
        if str(stas[ii])=='AID':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+2.838e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='ATD':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+2.840e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='BILL':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+4.305e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='BORG':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+3.17e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='CHTO':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+3.627e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='EFI':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+2.75e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='GUMO' and str(chans[ii])=='LHZ' and str(locs[ii])=='10':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+1.798e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='HRV':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+4.58e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='MIDW':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            #w(ii) = extract(w(ii),'TIME',stime+0.75e4/spdy,get(w(ii),'end'))
            w[ii].trim(stime+0.75e4,tr.stats.endtime)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='PFO':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+3.3205e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='SCZ':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+1.93e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='TAM':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+5.58e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='TIXI':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+6.615e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        if str(stas[ii])=='WUS':
            print('special cutting of the end of the record');
            stime = w[ii].stats.starttime
            w[ii].trim(stime,stime+5.82e5)
            if ifigure==1:
                figname=plt.figure()
                plt.plot(w[ii].data,'r') 
                plt.title(str(stag)+ ' special cut')
        
        
        # skip these
        if str(stas[ii])=='ADK' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='DAV' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='DGAR' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='GRFO' and str(chans[ii])=='LHZ':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='MBWA' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='RAO' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='POHA' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='WAKE' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='XMAS' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='HOPE' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='LCO' and str(chans[ii])=='LHZ':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='PTCN' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='PEL' and str(chans[ii])=='LHZ':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='RSSD' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='SHEL' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='PALK' and str(chans[ii])=='LHZ' and str(locs[ii])=='00':
            print('SKIPPING'); scut.append(ii)
        
        
        # skip these (multiple records per station)
        
        if str(stas[ii])=='TRISIU':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='TRIS':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='FUNA':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='OTAV':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='PMG':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='SAML':
            print('SKIPPING'); scut.append(ii)
        
        if str(stas[ii])=='SDV':
            print('SKIPPING'); scut.append(ii)
        

        # check for gaps and fill with mean data
        
        wd = w[ii].data
        if len(wd) !=0:
            values = np.array(wd)
            searchval = 0
            icut = np.where(values == searchval)[0]
            mnst=stats.mean(wd)
            if len(icut) != 0:
                if len(icut)==1:
                    print('filling the one-point gap');
                    w[ii].data = np.where(w[ii].data ==0.0000, mnst, w[ii].data)
                else:
                    fgap = len(icut)/len(wd);
                    if fgap < FTHRESH:
                        stgap = ('filling %.3f (< %.3f) gaps with meanAll' % (fgap,FTHRESH))
                        print(stgap)
                        w[ii].data = np.where(w[ii].data ==0.0000, mnst, w[ii].data)
                        if ifigure==1:
                            figname=plt.figure()
                            plt.plot(w[ii].data,'r') 
                            plt.title(str(stit)+ ': '+ str(stgap))
          
                    else:
                        print(fgap)
                        w[ii].data = np.where(w[ii].data ==0.0000, np.nan, w[ii].data)
                        print('FIX THIS');
                        print(['skipping '+ str(stag)+ ' due to gaps']);
                        #icut[ii] = [icut ii];
                        scut.append(ii)
        else:
            print('no data for station',stas[ii])
            scut.append(ii)
            # even number seems to make the FFT work easier -- cut a point if needed
            #wd = w[ii].data;
            #if mod(length(wd),2)==1
             #   print('cutting the last point to make it an even number (FFT)');
              #  w(ii) = extract(w(ii),'INDEX',1,length(wd)-1);
           

        ii+=1
    
    # write to file
    
    if savefile==True:
        filename = 'sumatra_modes.txt'
        print('writing %i points to file %s' % (nw,filename))
        with open(filename, "w") as file1:
            pdfp = -1
            for jj in range(nw):
                stag = [str(stas[jj])+ '_' +str(chans[jj])+ '_' +str(nets[jj])]
                if jj in scut:
                    pkeep=0
                else:
                    pkeep=1
                pdfp = pdfp + pkeep
                # Writing data to a file
                file1.write('%3i %3i %7s %7s %4s %16s %4i\n' %(jj,pdfp*pkeep,str(stas[jj]+str(locs[jj])),
                            str(chans[jj]),str(nets[jj]),str(stag),pkeep))
        
        file1.close()
    
    
    # compute fft for each time series that is NOT cut
    directory1 = "./datawf/sumatra_fft"
    if not path.exists(directory1):  # If data directory doesn't exist, it will create one
        os.makedirs(directory1)
    all_fft_f=[]
    all_fft_amp=[]
    all_fft_phs=[]
    ii=imin
    while ii < imax:
        if ii not in scut:
            stag = str(stas[ii])+ '_'+ str(chans[ii])+ '_' +str(nets[ii])
            dur = len(w[ii])/spdy;
            stdur = 'duration = %.2f days'% (dur)
            stit = str(stag)+ ', '+ str(stdur)
            print('%i/%i %s' % (ii,nw,stag))
            dur_s = 10*24*60*60    # convert 10 days into s
            
            tr = w[ii].copy()
            mnst=stats.mean(tr.data)
            t = tr.stats.starttime
           
            # Make sure they have the same no. of sample
            tr.trim(t,t+dur_s, pad=True, fill_value=mnst)
            
            w_detrend=detrend(tr.data,'constant')
            w_demean=w_detrend-mnst
            taper_percentage = 1
            npts = len(w_demean)              # number of samples
            df = tr.stats.sampling_rate       # sampling rate
            nsec = npts/df                        # sampling time
            fNy = df / 2.0                        # Nyquist frequency
            taper = cosine_taper(npts,taper_percentage)
            w_taper = w_demean * taper
            fft_amp, fft_phase, f = wf_fft(w_taper,fNy) # amplitude, and phase, frequencies of fft
            
            all_fft_f.append(f)
            all_fft_amp.append(fft_amp)
            all_fft_phs.append(fft_phase)
            
            # write and save the stream into ".sac" format
            tr_fft = Trace(np.array(fft_amp))
            tr_fft.write(path.join(directory1, stas[ii]+locs[ii])+"amps", format = 'SAC')  
        #print (tr_fft)
        ii+=1
    # save amplitudes and frequencies to npy dat file
    np.save(path.join(directory1, 'all_fft_freq'), all_fft_f)
    np.save(path.join(directory1, 'all_fft_phase'), all_fft_phs)            
    
    wnew=fft_amp
    return wnew 