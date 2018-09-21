# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 09:54:53 2018

@author: Alyssa
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

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
    
    # Now read in the actual noise data
    noisefile = folder + 'noiseSet0000'+ PnFn +'.bin'
    streamdata = np.fromfile(noisefile,'>f8') #first half of points are I, second half are Q
    split = int(len(streamdata)/2) # find the index of the halfway point
    streamrawI = streamdata[0:split]
    streamrawQ = streamdata[split:]
    d[Pn][Fn]['stream']['raw S21'] = u.V*streamrawI+1j*u.V*streamrawQ

    
def importmixerfolder(d,T_stage,T_BB,cool,folder='',docal=True,Qignore=10**3,poly_order=5):
    atten_list = np.loadtxt(folder+'attenuations.txt',delimiter=',')
    freq_list = np.loadtxt(folder+'initial_f0.txt',delimiter=',')
    
    for p in np.arange(0,len(atten_list)):
        for f in np.arange(0,len(freq_list)):
            importmixerdata(d,T_stage,T_BB,cool,folder=folder,Pn=p,Fn=f)
            if docal==True:
                gaincal(d[p][f],Qignore,poly_order)


def gaincal(resdict,Qignore=10**3,poly_order=5):
    
    try: f00 = resdict['stream']['freq'] # If there's noise data, use the frequency at which it was streamed as the center frequency
    except KeyError: f00 = resdict['initial f0'] # If there's no noise data, just use the initial input frequency as the center frequency
    df = f00/Qignore
    
    f = resdict['gain']['freqs']
    rawS21 = resdict['gain']['raw S21']

    conds = (f < f00-df)|(f > f00+df) #for the fits, ignore the area within df of the resonance
    
    #Fit & divide out any gain variations across the band, ignoring the resonance itself
    amp_poly = np.polyfit(f[conds],np.abs(rawS21[conds]),poly_order)
    amp_cal_gain = np.polyval(amp_poly,f)
    
    #Fit & divide out the cable delay
    phi = np.arctan2(np.imag(rawS21),np.real(rawS21))
    
    cable_lin = np.polyfit(f[conds],phi[conds],1)
    phi_cal_gain = np.polyval(cable_lin,f) # this is tau*f + phi_0
    
    resdict['gain']['cal S21'] = rawS21*np.exp(-1j*phi_cal_gain)/amp_cal_gain
    
    resdict['gain']['phi cal gain'] = phi_cal_gain
    resdict['gain']['amp cal gain'] = amp_cal_gain
    
    # Now apply the gain calibrations to the fine + med + rough scans
    types = ['fine','med','rough']
    for t in types:
        amp_cal = np.polyval(amp_poly,resdict[t]['freqs']) # plug in the amp gain polynomial fit for the relevant frequency range
        phi_cal = np.polyval(cable_lin,resdict[t]['freqs']) # plug in the cable delay [linear] fit for the relevant frequency range
        resdict[t]['cal S21'] = resdict[t]['raw S21']*np.exp(-1j*phi_cal)/amp_cal # make the corrections for this cal scan type
    
    # Apply the gain calibration to the noise timestream
    amp_cal_noise = np.polyval(amp_poly,resdict['stream']['freq']) # plug in the amp gain polynomial fit for the relevant frequency point
    phi_cal_noise = np.polyval(cable_lin,resdict['stream']['freq']) # plug in the cable delay [linear] fit for the relevant frequency point
    resdict['stream']['cal S21'] = resdict['stream']['raw S21']*np.exp(-1j*phi_cal_noise)/amp_cal_noise # correct the noise stream
        
       
#%%
# From https://gist.github.com/lorenzoriano/6799568
def calc_R(x,y, xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f(c, x, y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()

from scipy import optimize
def leastsq_circle(x,y):
    # coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, args=(x,y))
    xc, yc = center
    Ri       = calc_R(x, y, *center)
    R        = Ri.mean()
    residu   = np.sum((Ri - R)**2)
    return xc, yc, R, residu

#%%

def freq_phase_func(freq,f0,Qr,th00):
    th = -th00 + 2*np.arctan(2*Qr*(1-freq/f0))
    return th

from scipy.optimize import curve_fit
def fit_freq_phase(freqs,S21,f00):
    phases = np.unwrap(np.angle(S21))
    
    popt,pcov = curve_fit(freq_phase_func,freqs,phases,p0=(f00,1e5,0))
    return popt,pcov

#%%
from scipy.interpolate import CubicSpline
def streamcal(resdict,phase_poly_order):
    calI = np.real(resdict['fine']['cal S21']).value
    calQ = np.imag(resdict['fine']['cal S21']).value
    
    # fit the IQ loop to a circle
    xc,yc,R,resid = leastsq_circle(calI,calQ)
    
    # shift the center of the circle to the origin
    cor_calI = calI - xc
    cor_calQ = calQ - yc
    
    # find the calibration point corresponding to the streaming freq point
    freqs = resdict['fine']['freqs']
    f00 = resdict['stream']['freq']
    
    ind = np.where(freqs==f00)[0] 
    if ind.size==0: ind = np.max(np.where(freqs<f00)) # not sure if the mixer always chooses an exact point or not
    
    # find the angle between the +I axis and the streaming freq point
    i0 = cor_calI[ind]
    q0 = cor_calQ[ind]
    th0 = np.unwrap(np.angle(i0+1j*q0))
    
    # rotate the IQ loop so that the streaming freq point lies on the +I axis
    rot_cor_calS21 = (cor_calI+1j*cor_calQ)*np.exp(-1j*th0)
    phase = np.unwrap(np.angle(rot_cor_calS21))
    
    
    # apply the shift and rotation to the noise stream
    stream_calI = np.real(resdict['stream']['cal S21']).value
    stream_calQ = np.imag(resdict['stream']['cal S21']).value
    
    stream_cor_calI = stream_calI - xc
    stream_cor_calQ = stream_calQ - yc
    
    stream_rot_cor_calS21 = (stream_cor_calI+1j*stream_cor_calQ)*np.exp(-1j*th0)
    
    # find the phase range over which we have streaming data
    stream_rot_cor_phase = np.unwrap(np.angle(stream_rot_cor_calS21))
    
    popt,pcov = fit_freq_phase(freqs.value,rot_cor_calS21,f00.value)
    f0,Qr,th0 = popt
    
    Qc = Qr*(np.sqrt(xc**2+yc**2)+R)/(2*R)
    
    phasesort = np.argsort(phase)
    freq_phase_interp = CubicSpline(phase[phasesort],freqs[phasesort])
    stream_freqs = freq_phase_interp(stream_rot_cor_phase)
    
    x_noise = (stream_freqs-f0)/f0
    sxx_psd,sxx_freqs = plt.psd(x_noise,Fs=20000,NFFT=2**14)
    
                    

#%%
    
testdict = {}
folder = '20180607/Noise01/'
T_stage = 0.215*u.K
T_BB = 5.61*u.K
cool = 'CD012'
importmixerfolder(testdict,T_stage,T_BB,cool,folder,docal=True,Qignore=10**3,poly_order=5) 

t1 = testdict[0][0]
f1 = t1['gain']['freqs']
rawI1 = np.real(t1['gain']['raw S21'])
rawQ1 = np.imag(t1['gain']['raw S21'])
calI1 = np.real(t1['gain']['cal S21'])
calQ1 = np.imag(t1['gain']['cal S21'])

plt.plot(rawI1,rawQ1,'r.')
plt.plot(calI1,calQ1,'b.')
    
rawI2 = np.real(t1['fine']['raw S21'])
rawQ2 = np.imag(t1['fine']['raw S21'])
calI2 = np.real(t1['fine']['cal S21'])
calQ2 = np.imag(t1['fine']['cal S21'])

plt.plot(rawI2,rawQ2,'m.')
plt.plot(calI2,calQ2,'c.')
plt.plot(1,0,'kx')

xc,yc,R,resid = leastsq_circle(calI2.value,calQ2.value)

from matplotlib.patches import Circle

fig, ax = plt.subplots()
ax.add_artist(Circle((xc,yc),R,fill=False))

    
