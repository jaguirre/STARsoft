# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 09:54:53 2018

@author: Alyssa
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from AnalyzeResonatorData import *

import warnings
warnings.simplefilter('ignore',np.RankWarning)

def importmixerdata(d,T_stage,T_BB,cool,folder='',Pn=0,Fn=0):
    Pn_label = str(Pn).zfill(2)
    Fn_label = str(Fn).zfill(2)
    PnFn = 'Pn'+Pn_label+'Fn'+Fn_label

    # should add something here to indicate whether d[Pn][Fn] already has data in it
    try: d[Pn][Fn] = {}
    except KeyError: 
        d[Pn] = {}
        d[Pn][Fn] = {}
    
    atten_list = np.loadtxt(folder+'attenuations.txt',delimiter=',')
    d[Pn][Fn]['LB_atten'] = u.dB*atten_list[Pn]
    
    freq_list = np.loadtxt(folder+'initial_f0.txt',delimiter=',')
    d[Pn][Fn]['initial f0'] = freq_list[Fn]
    
    d[Pn][Fn]['T_BB'] = T_BB
    d[Pn][Fn]['T_stage'] = T_stage
    d[Pn][Fn]['folder'] = folder
    d[Pn][Fn]['CD'] = cool
    
    # read in the calibration data
    types = ['fine','med','gain','rough']
    for t in types:
        sweepfile = folder + t + 'SweepSet0000' + PnFn + '.txt'
        try: freqs,rawI,rawQ = np.loadtxt(sweepfile,delimiter=',',skiprows=3,unpack=True)
        except OSError: print('No cal data found for ' + t + '--' + PnFn + '; moving on')
        else:
            d[Pn][Fn][t] = {}
            d[Pn][Fn][t]['freqs'] = u.Hz*freqs
            d[Pn][Fn][t]['raw S21'] = u.V*rawI+1j*u.V*rawQ

    # grab the frequency at which the noise was streamed
    noisefreqfile = folder + 'noiseFrequencySet0000' + PnFn + '.txt'
    t1 = open(noisefreqfile)
    noisefreq = u.Hz*float(t1.readlines()[2][21:-1]) # hard-coding in, merp. there must be a better way...
    t1.close()
    d[Pn][Fn]['stream'] = {}
    d[Pn][Fn]['stream']['freq'] = noisefreq
    
    # also grab the data acquisition rate of the mixer system
    logfile = folder + 'log.txt'
    t2 = open(logfile)
    streamrate = u.Hz*float(t2.readlines()[4][17:-1]) # hard-coding in, merp. there must be a better way...
    t2.close()
    d[Pn][Fn]['stream']['streamrate'] = streamrate
    
    # Now read in the actual noise data
    noisefile = folder + 'noiseSet0000'+ PnFn +'.bin'
    streamdata = np.fromfile(noisefile,'>f8') #first half of points are I, second half are Q
    split = int(len(streamdata)/2) # find the index of the halfway point
    streamrawI = streamdata[0:split]
    streamrawQ = streamdata[split:]
    d[Pn][Fn]['stream']['raw S21'] = u.V*streamrawI+1j*u.V*streamrawQ

    
def importmixerfolder(d,T_stage,T_BB,cool,folder='',docal=True,Qignore=10**3,poly_order=5,doPSD=True):
    atten_list = np.loadtxt(folder+'attenuations.txt',delimiter=',')
    freq_list = np.loadtxt(folder+'initial_f0.txt',delimiter=',')
    
    for p in np.arange(0,len(atten_list)):
        for f in np.arange(0,len(freq_list)):
            importmixerdata(d,T_stage,T_BB,cool,folder=folder,Pn=p,Fn=f)
            if docal==True:
                gaincal(d[p][f],Qignore,poly_order)
                streamcal(d[p][f])
            if doPSD==True:
                streampsd(d[p][f])
                
                


#%%
testdict = {}
folder = '20180607/Noise01/'
T_stage = 0.215*u.K
T_BB = 5.61*u.K
cool = 'CD012'
importmixerfolder(testdict,T_stage,T_BB,cool,folder,docal=True,Qignore=10**3,poly_order=5) 

#%%
resdict = testdict[2][2]
plt.figure(-6)
comb_freqs = np.concatenate(((resdict['med']['freqs'].to(u.Hz)).value,(resdict['fine']['freqs'].to(u.Hz)).value))
inds = np.argsort(comb_freqs)
comb_s21 = np.concatenate((resdict['med']['cal S21'],resdict['fine']['cal S21']))
freqs = u.Hz*comb_freqs[inds]
S21 = comb_s21[inds]
plt.plot(freqs,np.abs(S21),'ko',label='med and fine data')
Qrt = resdict['Qr_calc']
Qct = resdict['Qc_calc']
f0t = resdict['f0_calc']
S21t = S21_linear(freqs.value,f0t,Qrt,Qct)
plt.plot(freqs,np.abs(S21t),'m--',label='cal simple model')

it = S21.real
qt = S21.imag
S21tofit = np.concatenate((it,qt))

p0t = (f0t,Qrt,Qct,0,0)
boundst = ([resdict['fine']['freqs'].min().value,1e3,1e3,-np.inf,-1],[resdict['fine']['freqs'].max().value,1e6,1e8,np.inf,1])
res_popt,res_pcov = curve_fit(fit_S21_nonlinear,freqs.value,S21tofit,p0=p0t,bounds=boundst)
plt.plot(freqs,np.abs(S21_nonlinear(freqs.value,*res_popt)),'c-',label='nonlinear S21 fit')
plt.xlabel('freq (Hz)')
plt.ylabel('|S21|')
plt.legend(loc='lower left')
#%%
#chi = 0
#plt.figure('test grr')
#plt.plot(freqs,np.abs(S21_nonlinear(freqs,f0t,Qrt,Qct,chi,0)),'b.')
##%%
#
#for sweep,scan in testdict:
#    resdict = testdict[sweep][scan]
#    comb_freqs = np.concatenate(((resdict['med']['freqs'].to(u.Hz)).value,(resdict['fine']['freqs'].to(u.Hz)).value))
#    inds = np.argsort(comb_freqs)
#    comb_s21 = np.concatenate((resdict['med']['cal S21'],resdict['fine']['cal S21']))
#    freqs = u.Hz*comb_freqs[inds]
#    S21 = comb_s21[inds]
#    plt.plot(freqs,np.abs(S21),'ko',label='med and fine data')
#    Qrt = resdict['Qr_calc']
#    Qct = resdict['Qc_calc']
#    f0t = resdict['f0_calc']
#    S21t = S21_linear(freqs.value,f0t,Qrt,Qct)
#    plt.plot(freqs,np.abs(S21t),'m--',label='cal simple model')
#    
#    it = S21.real
#    qt = S21.imag
#    S21tofit = np.concatenate((it,qt))
#    
#    p0t = (f0t,Qrt,Qct,0,0)
#    boundst = ([resdict['fine']['freqs'].min().value,1e3,1e3,-np.inf,-1],[resdict['fine']['freqs'].max().value,1e6,1e8,np.inf,1])
#    res_popt,res_pcov = curve_fit(fit_S21_nonlinear,freqs.value,S21tofit,p0=p0t,bounds=boundst)
#    plt.plot(freqs,np.abs(S21_nonlinear(freqs.value,*res_popt)),'c-',label='nonlinear S21 fit')
#    plt.xlabel('freq (Hz)')
#    plt.ylabel('|S21|')
#    plt.legend(loc='lower left')
#    
#    a = res_popt
#    plt.title('fit vals: ' + 'f0 = ' + str(a[0]) + ', Qr = ' + str(a[1]) + ', Qc = ' + str(a[2]) + ', chi = ' + str(a[3]) + ', anl = ' + str(a[4]),fontsize=8)



#%%
#        
#t1 = testdict[0][0]
#f1 = t1['gain']['freqs']
#rawI1 = np.real(t1['gain']['raw S21'])
#rawQ1 = np.imag(t1['gain']['raw S21'])
#calI1 = np.real(t1['gain']['cal S21'])
#calQ1 = np.imag(t1['gain']['cal S21'])
#
#plt.plot(rawI1,rawQ1,'r.')
#plt.plot(calI1,calQ1,'b.')
#    
#rawI2 = np.real(t1['fine']['raw S21'])
#rawQ2 = np.imag(t1['fine']['raw S21'])
#calI2 = np.real(t1['fine']['cal S21'])
#calQ2 = np.imag(t1['fine']['cal S21'])
#
#plt.plot(rawI2,rawQ2,'m.')
#plt.plot(calI2,calQ2,'c.')
#plt.plot(1,0,'kx')
#
#xc,yc,R,resid = leastsq_circle(calI2.value,calQ2.value)
#
#from matplotlib.patches import Circle
#
#fig, ax = plt.subplots()
#ax.add_artist(Circle((xc,yc),R,fill=False))
#
#    
