# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 10:15:03 2018

@author: Alyssa
"""
import numpy as np
from astropy import units as u

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
           