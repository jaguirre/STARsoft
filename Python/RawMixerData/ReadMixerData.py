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
from astropy.stats import sigma_clip
import scipy.signal as sig
def deglitch(S21,clip_sigma=5,clip_iters=5):
    I_orig = S21.real
    Q_orig = S21.imag
    
    I_sc = sigma_clip(I_orig,sigma=clip_sigma,iters=clip_iters)
    Q_sc = sigma_clip(Q_orig,sigma=clip_sigma,iters=clip_iters)
    
    good_data = ~np.logical_or(I_sc.mask,Q_sc.mask)
    S21_deglitched = I_orig[good_data] + 1j*Q_orig[good_data]
    
    return S21_deglitched

#%%
from scipy.interpolate import CubicSpline
def streamcal(resdict):
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
    ph0 = phase[ind]
    if np.abs(ph0) > np.pi: phase = phase-ph0
    
    
    # apply the shift and rotation to the noise stream
    stream_calI = np.real(resdict['stream']['cal S21']).value
    stream_calQ = np.imag(resdict['stream']['cal S21']).value
    
    stream_cor_calI = stream_calI - xc
    stream_cor_calQ = stream_calQ - yc
    
    stream_rot_cor_calS21_raw = (stream_cor_calI+1j*stream_cor_calQ)*np.exp(-1j*th0)
    resdict['stream']['rot cor cal S21 raw'] = stream_rot_cor_calS21_raw
    
    stream_rot_cor_calS21_degl = deglitch(stream_rot_cor_calS21_raw)
    resdict['stream']['rot cor cal S21 deglitched'] = stream_rot_cor_calS21_degl
    
    # find the phases of the streaming data
    stream_rot_cor_phase = np.unwrap(np.angle(stream_rot_cor_calS21_degl))
    
    # fit the phase vs freq curve to a function
    popt,pcov = fit_freq_phase(freqs.value,rot_cor_calS21,f00.value)
    # the fit values are f0 and Qr in a simplified resonator model
    f0,Qr,th0 = popt 
    
    # calculate what Qc would be in the simplified resonator model
    Qc = Qr*(np.sqrt(xc**2+yc**2)+R)/(2*R)
    
    # save the fit values as initial guesses for a more complex model
    resdict['f0_calc'] = f0
    resdict['Qr_calc'] = Qr
    resdict['Qc_calc'] = Qc
    
    # ensure that the phases are in order from least to greatest (could be problematic if there's a glitch??)
    phasesort = np.argsort(phase)
    
    # interpolate to get a function to go from phase -> freq
    freq_phase_interp = CubicSpline(phase[phasesort],freqs[phasesort])
    
    # plug the streaming points into the conversion function to get a frequency timestream
    stream_freq_noise = freq_phase_interp(stream_rot_cor_phase)
    
    # use the f0 value to convert to fractional frequency noise
    x_noise = (stream_freq_noise-f0)/f0
    resdict['stream']['x noise'] = x_noise

#%%    
def whitenoise(noise_spectrum,white_freq_range):
    freqs = noise_spectrum[0]
    psd = noise_spectrum[1]
    
    start_freq = white_freq_range[0]
    start_ind = np.max(np.where(freqs<start_freq)) + 1

    stop_freq = white_freq_range[1]
    stop_ind = np.min(np.where(freqs>stop_freq))
    
    white_level = np.average(psd[start_ind:stop_ind])
    return white_level

#%% 
from matplotlib.mlab import psd as psdnoplot
invHz = np.power(u.Hz,-1)

def streampsd(resdict,white_freq_range=[30,100]*u.Hz):
    stream_rot_cor_calS21 = resdict['stream']['rot cor cal S21 deglitched']
    stream_x_noise = resdict['stream']['x noise']
    streamrate = ((resdict['stream']['streamrate']).to(u.Hz)).value
    
    
    # calculate the parallel and perpendicular contributions to the noise
    perp_psd,perp_freqs = psdnoplot(stream_rot_cor_calS21.real,Fs=streamrate,NFFT=2**14,scale_by_freq=True)
    resdict['stream']['perp'] = [perp_freqs*u.Hz,perp_psd*invHz]

    par_psd,par_freqs = psdnoplot(stream_rot_cor_calS21.imag,Fs=streamrate,NFFT=2**14,scale_by_freq=True)
    resdict['stream']['par'] = [par_freqs*u.Hz,par_psd*invHz]
    
    # approximate the amplifier/electronic noise contribution to the total measured nosie
    rho = (par_psd-perp_psd)/par_psd
    
    # calculate the sxx psd 
    sxx_raw_psd,sxx_freqs = psdnoplot(stream_x_noise,Fs=streamrate,NFFT=2**14,scale_by_freq=True)
    resdict['stream']['raw Sxx'] = [sxx_freqs*u.Hz,sxx_raw_psd*invHz]
    
    # correct for the amplifier contribution
    resdict['stream']['amp sub Sxx'] = [sxx_freqs[1:]*u.Hz,rho[1:]*sxx_raw_psd[1:]*invHz]

    # calculate the white noise levels
    resdict['stream']['raw Sxx white'] = whitenoise(resdict['stream']['raw Sxx'],white_freq_range)
    resdict['stream']['amp sub Sxx white'] = whitenoise(resdict['stream']['amp sub Sxx'],white_freq_range)
    
#%%
testdict = {}
folder = '20180607/Noise01/'
T_stage = 0.215*u.K
T_BB = 5.61*u.K
cool = 'CD012'
importmixerfolder(testdict,T_stage,T_BB,cool,folder,docal=True,Qignore=10**3,poly_order=5) 

#%%
    
plt.figure(-1)
colors = ['r','g','b','c']
for pn in np.arange(0,9):
#    print('pn = ' + str(pn))
    for fn in np.arange(0,4):
        plt.plot(testdict[pn][fn]['LB_atten'],testdict[pn][fn]['stream']['raw Sxx white'],'+',color=colors[fn])
#        print('fn = ' + str(fn))
#        print(testdict[pn][fn]['stream']['raw Sxx white'])
        
#%%
plt.figure(-2)
comb_freqs = np.concatenate((resdict['med']['freqs'],resdict['fine']['freqs']))
inds = np.argsort(comb_freqs)
comb_s21 = np.concatenate((resdict['med']['cal S21'],resdict['fine']['cal S21']))
freqs = comb_freqs[inds]
S21 = comb_s21[inds]
plt.plot(freqs,np.abs(S21),'ko')
Qrt = resdict['Qr_calc']
Qct = resdict['Qc_calc']
f0t = resdict['f0_calc']
S21t = S21_linear(freqs.value,f0t,Qrt,Qct)
plt.plot(freqs,np.abs(S21t),'m--')

it = S21.real
qt = S21.imag
S21t = np.concatenate((it,qt))

p0t = (f0t,Qrt,Qct,0)
boundst = ([resdict['fine']['freqs'].min().value,1e3,1e3,-1],[resdict['fine']['freqs'].max().value,1e6,1e6,1])
res_popt,res_pcov = curve_fit(fit_S21_nonlinear,freqs.value,S21t,p0=p0t,bounds=boundst)
plt.plot(freqs,np.abs(S21_nonlinear(freqs,*res_popt)),'c-')
#%%
def S21_linear(fr,f0,Qr,Qc):
    x = (fr-f0)/f0
    S21 = 1-(Qr/Qc)*np.power(1+2*1j*Qr*x,-1)
    return S21

#def S21_nonlinear(fr,f0,Qr,Qc,chi,a_nl):
#    y0 = (fr-f0)/f0
#    
#    # solution to y=y0+a_nl/(1+4*y**2) as stated in McCarrick 2014 appendix (original reference is Swenson 2013):
#    k2 = np.power(((y0**3/27 + y0/12 + a_nl/8)**2 - (y0**2/9-1/12)**3),1/2)
#    
#    k1 = np.power((a_nl/8 + y0/12 + k2 + y0**3/27),1/3)
#    
#    y = y0/3 + ((y0**2/9 - 1/12)/k1) + k1
#    
#    x = y/Qr 
#    
#    S21 = 1-(Qr/Qc)*(1+1j*chi)*np.power(1+2*1j*Qr*x,-1)
#    return S21

def S21_nonlinear(fr,f0,Qr,Qe,a_nl):
    y0 = (fr-f0)/f0
    
    # solution to y=y0+a_nl/(1+4*y**2) as stated in McCarrick 2014 appendix (original reference is Swenson 2013):
    k2 = np.power(((y0**3/27 + y0/12 + a_nl/8)**2 - (y0**2/9-1/12)**3),1/2)
    
    k1 = np.power((a_nl/8 + y0/12 + k2 + y0**3/27),1/3)
    
    y = y0/3 + ((y0**2/9 - 1/12)/k1) + k1
    
    x = y/Qr 
    
    S21 = 1-(Qr/Qe)*np.power(1+2*1j*Qr*x,-1)
    return S21


def I_nonlinear(fr,f0,Qr,Qc,a_nl):
    S21 = S21_nonlinear(fr,f0,Qr,Qc,a_nl)
    I = S21.real
    return I

def Q_nonlinear(fr,f0,Qr,Qc,a_nl):
    S21 = S21_nonlinear(fr,f0,Qr,Qc,a_nl)
    Q = S21.imag
    return Q

def fit_S21_nonlinear(fr,f0,Qr,Qc,a_nl):
    I = I_nonlinear(fr,f0,Qr,Qc,a_nl)
    Q = Q_nonlinear(fr,f0,Qr,Qc,a_nl)
    
    glom = np.concatenate((I,Q))
    return glom
#%%
        
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

    
