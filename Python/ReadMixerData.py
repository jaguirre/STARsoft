# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 09:54:53 2018

@author: Alyssa
"""

import numpy as np
import matplotlib as mpl
mpl.rcParams['grid.linestyle'] = ':'
import matplotlib.pyplot as plt
from astropy import units as u
from AnalyzeResonatorData import *

import warnings
warnings.simplefilter('ignore',np.RankWarning)

def importmixerdata(d,T_stage,T_BB,cool,datafolder='',outfolder='',datescan='',Pn=0,Fn=0,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5):
    Pn_label = str(Pn).zfill(2)
    Fn_label = str(Fn).zfill(2)
    PnFn = 'Pn'+Pn_label+'Fn'+Fn_label
    print('Importing mixer data ' + cool + ' :' + datescan + PnFn)
    # should add something here to indicate whether d[Pn][Fn] already has data in it
    try: d[Pn][Fn] = {}
    except KeyError: 
        d[Pn] = {}
        d[Pn][Fn] = {}
    
    atten_list = np.loadtxt(datafolder+'attenuations.txt',delimiter=',')
    d[Pn][Fn]['LB_atten'] = u.dB*atten_list[Pn]
    
    freq_list = np.loadtxt(datafolder+'initial_f0.txt',delimiter=',')
    try: d[Pn][Fn]['initial f0'] = freq_list[Fn]
    except IndexError: d[Pn][Fn]['initial f0'] = freq_list
    
    d[Pn][Fn]['T_BB'] = T_BB
    d[Pn][Fn]['T_stage'] = T_stage
    d[Pn][Fn]['folder'] = datafolder
    d[Pn][Fn]['CD'] = cool
    d['CD'] = cool
    
    # read in the calibration data
    types = ['fine','med','gain','rough']
    for t in types:
        sweepfile = datafolder + t + 'SweepSet0000' + PnFn + '.txt'
        try: freqs,rawI,rawQ = np.loadtxt(sweepfile,delimiter=',',skiprows=3,unpack=True)
        except OSError: print('No cal data found for ' + t + '--' + PnFn + '; moving on')
        else:
            d[Pn][Fn][t] = {}
            d[Pn][Fn][t]['freqs'] = u.Hz*freqs
            d[Pn][Fn][t]['raw S21'] = u.V*rawI+1j*u.V*rawQ

    # grab the frequency at which the noise was streamed
    noisefreqfile = datafolder + 'noiseFrequencySet0000' + PnFn + '.txt'
    t1 = open(noisefreqfile)
    noisefreq = u.Hz*float(t1.readlines()[2][21:-1]) # hard-coding in, merp. there must be a better way...
    t1.close()
    d[Pn][Fn]['stream'] = {}
    d[Pn][Fn]['stream']['freq'] = noisefreq
    
    # also grab the data acquisition rate of the mixer system
    logfile = datafolder + 'log.txt'
    t2 = open(logfile)
    streamrate = u.Hz*float(t2.readlines()[4][17:-1]) # hard-coding in, merp. there must be a better way...
    t2.close()
    d[Pn][Fn]['stream']['streamrate'] = streamrate
    
    # Now read in the actual noise data
    noisefile = datafolder + 'noiseSet0000'+ PnFn +'.bin'
    streamdata = np.fromfile(noisefile,'>f8') #first half of points are I, second half are Q
    split = int(len(streamdata)/2) # find the index of the halfway point
    streamrawI = streamdata[0:split]
    streamrawQ = streamdata[split:]
    d[Pn][Fn]['stream']['raw S21'] = u.V*streamrawI+1j*u.V*streamrawQ

    ######## Saving data
    # if the folder to save data in doesn't exist, make it
    
    if os.path.isdir(outfolder) == False: os.mkdir(outfolder)

    if docal==True:
        gaincal(d[Pn][Fn],Qignore,poly_order)
        streamcal(d[Pn][Fn])
        fit_resonator_S21(d[Pn][Fn])
    if doPSD==True:
        streampsd(d[Pn][Fn])
    if doplots==True:
        plt.close('all')
        
        Pn_label = str(Pn).zfill(2)
        Fn_label = str(Fn).zfill(2)
        PnFn = 'Pn'+Pn_label+'Fn'+Fn_label

        resdict = d[Pn][Fn]
        # plot freq vs phase conversion step
        fig1name = outfolder+'freq-phase-plot_'+cool+'_'+datescan+PnFn+'.png'
        fig1 = plt.figure(fig1name)
        p1 = fig1.add_subplot(111)
        p1.plot(*resdict['freq-phase plot']['fine data (rot cor cal)'],'o',color='k',label='fine data (rot cor cal)')
        p1.plot(*resdict['freq-phase plot']['fine fit to function'],'-',linewidth=2,color='r',label='fine fit to function')
        p1.plot(*resdict['freq-phase plot']['x noise'],'.',color='b',label='frequency noise')
#                p1.axhline(y=resdict['f0_calc'],linestyle='--',color='g',alpha=0.5,label='f0_calc')
        p1.plot(*resdict['freq-phase plot']['streaming freq'],'*',color='magenta',markersize=9,label='streaming freq')
        p1.legend(loc='lower left',frameon=True)
        p1.set_xlabel(r'$\delta$f/f',fontsize=12)
        p1.set_ylabel(r'$\theta$ (radians)',fontsize=12)
        p1.set_title(cool+'_'+datescan+PnFn)
        p1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        fig1.savefig(fig1name)
        plt.close(fig1name)
        
        # plot resonator line profile fit step
        fig2name = outfolder+'mag-S21-plot_'+cool+'_'+datescan+PnFn+'.png'
        fig2, (p2,p2b) = plt.subplots(2,1,gridspec_kw={'height_ratios':[3, 1]}, sharex=True, num=fig2name)
#        fig2 = plt.figure(fig2name)
#        p2 = fig2.add_subplot(211)
        p2.plot(*resdict['mag-S21-plot']['cal fit data'],'o',color='k',label='fine & med cal data')
        p2.plot(*resdict['mag-S21-plot']['nonlinear func fit'],'-',linewidth=2,color='r',label='nonlinear function fit')
        p2.legend(loc='lower left',frameon=True)
        p2.set_ylabel(r'|$S_{21}$| (dB)',fontsize=12)
        p2.set_title(cool+'_'+datescan+PnFn)
        p2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        fits = r'$f_0$ = {:.2f}'.format(resdict['f0_fit'].to(u.MHz)) + '\n' + r'$Q_r$ = {:.1e}'.format(resdict['Qr_fit']) + '\n' + r'$Q_c$ = {:.1e}'.format(resdict['Qc_fit']) + '\n' + r'$a_{nl}$ = ' +'{:.2f}'.format(resdict['anl_fit'])
        p2.text(0.95, 0.01, fits,verticalalignment='bottom',horizontalalignment='right',transform=p2.transAxes,color='k',fontsize=10)
        
#        p2b = fig2.add_subplot(212,sharex=p2)
        p2b.axhline(y=0,color='k',linestyle=':',alpha=0.5)
        p2b.plot(resdict['mag-S21-plot']['cal fit data'][0],resdict['mag-S21-plot']['cal fit data'][1]-resdict['mag-S21-plot']['nonlinear func fit'][1],'-',color='darkred')
        p2b.set_xlabel(r'$\delta$f/f',fontsize=12)
        p2b.set_ylabel(r'resid(|$S_{21}$|)',fontsize=12)
        fig2.subplots_adjust(hspace=0)
        fig2.savefig(fig2name)
        plt.close(fig2name)
        
        # plot PSDs and Sxx
        fig3name = outfolder+'PSD-plot_'+cool+'_'+datescan+PnFn+'.png'
        fig3 = plt.figure(fig3name)
        p3top = fig3.add_subplot(211)
        p3top.plot(*resdict['stream']['perp'],'-',linewidth=2,color='r',label='Perpendicular')
        p3top.plot(*resdict['stream']['par'],'-',linewidth=2,color='m',label='Parallel')
        p3top.set_yscale('log')
        p3top.set_ylim([1e-9,1e-6])
        p3top.set_ylabel(r'$\partial S_{21}$ (Hz$^{-1}$)',fontsize=12)
        p3top.set_xscale('log')
        p3top.legend(loc='upper right',frameon=True,framealpha=1)
        p3top.tick_params(labelbottom=False)
        p3top.set_title(cool+'_'+datescan+PnFn)
        
        p3bot = fig3.add_subplot(212,sharex=p3top)
        p3bot.plot(*resdict['stream']['raw Sxx'],'-',linewidth=2,color='g',label='raw')
        p3bot.plot(*resdict['stream']['amp sub Sxx'],'-',linewidth=2,color='darkblue',label='amp. sub.')
        p3bot.plot(resdict['stream']['white_freq_range'],np.ones_like(resdict['stream']['white_freq_range'])*resdict['stream']['raw Sxx white'],'-',linewidth=1.5,color='limegreen',label='raw white')
        p3bot.plot(resdict['stream']['white_freq_range'],np.ones_like(resdict['stream']['white_freq_range'])*resdict['stream']['amp sub Sxx white'],'-',linewidth=1.5,color='dodgerblue',label='amp. sub. white')
        p3bot.set_yscale('log')
        p3bot.set_ylim([1e-18,1e-15])
        p3bot.set_ylabel(r'$S_{xx}$ (Hz$^{-1}$)',fontsize=12)
        p3bot.set_xscale('log')
        p3bot.set_xlabel('frequency (Hz)',fontsize=12)
        p3bot.legend(loc='upper right',frameon=True,framealpha=1)
        
        fig3.savefig(fig3name)
        plt.close(fig3name)
        
        
import os    
def importmixerfolder(d,T_stage,T_BB,cool,datafolder='',outfolder='',datescan='',docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5):
    atten_list = np.loadtxt(datafolder+'attenuations.txt',delimiter=',')
    freq_list = np.loadtxt(datafolder+'initial_f0.txt',delimiter=',')
    

    for p in np.arange(0,len(atten_list)):
        try:
            for f in np.arange(0,len(freq_list)):
                Pn_label = str(p).zfill(2)
                Fn_label = str(f).zfill(2)
                PnFn = 'Pn'+Pn_label+'Fn'+Fn_label
    
                if os.path.isfile(datafolder+'fineSweepSet0000'+PnFn+'.txt') and os.path.isfile(datafolder+'gainSweepSet0000'+PnFn+'.txt') and os.path.isfile(datafolder+'noiseSet0000'+PnFn+'.bin'): # for cases where the labview aborted in the middle of a run
                    importmixerdata(d,T_stage,T_BB,cool,datafolder,outfolder,datescan,p,f,docal,doPSD,doplots)

        except TypeError: #if there's only one frequency in the scan
            f = 0
            Pn_label = str(p).zfill(2)
            Fn_label = str(f).zfill(2)
            PnFn = 'Pn'+Pn_label+'Fn'+Fn_label

            if os.path.isfile(datafolder+'fineSweepSet0000'+PnFn+'.txt') and os.path.isfile(datafolder+'gainSweepSet0000'+PnFn+'.txt') and os.path.isfile(datafolder+'noiseSet0000'+PnFn+'.bin'): # for cases where the labview aborted in the middle of a run
                importmixerdata(d,T_stage,T_BB,cool,datafolder,outfolder,datescan,p,f,docal,doPSD,doplots)

#import deepdish as dd
from DictionaryToHDF5 import *
def savedatadict(testdict,cool,datescan,outfolder):
    if os.path.isdir(outfolder) == False: os.mkdir(outfolder)
    

    dictname = outfolder+cool+'_'+datescan+'_datadict.hdf5'
    dict_to_hdf5(testdict,dictname)
#    dd.io.save(dictname,testdict)
    


#%%
#testdict = {}
#datafolder = '20180607/Noise01/'
#T_stage = 0.215*u.K
#T_BB = 5.61*u.K
#cool = 'CD012'
#outfolder = 'C:/Users/Alyssa/Penn Google Drive/Penn & NSTRF/Caltech Devices/' + 'test/'
##importmixerdata(testdict,T_stage,T_BB,cool,datafolder,outfolder,Pn=1,Fn=1,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5)
#importmixerfolder(testdict,T_stage,T_BB,cool,datafolder,outfolder,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5) 
#
#savedatadict(testdict,cool,outfolder)
#
#
#


#%%
#resdict = testdict[2][2]
#plt.figure(-6)
#comb_freqs = np.concatenate(((resdict['med']['freqs'].to(u.Hz)).value,(resdict['fine']['freqs'].to(u.Hz)).value))
#inds = np.argsort(comb_freqs)
#comb_s21 = np.concatenate((resdict['med']['cal S21'],resdict['fine']['cal S21']))
#freqs = u.Hz*comb_freqs[inds]
#S21 = comb_s21[inds]
#plt.plot(freqs,np.abs(S21),'ko',label='med and fine data')
#Qrt = resdict['Qr_calc']
#Qct = resdict['Qc_calc']
#f0t = resdict['f0_calc']
#S21t = S21_linear(freqs.value,f0t,Qrt,Qct)
#plt.plot(freqs,np.abs(S21t),'m--',label='cal simple model')
#
#it = S21.real
#qt = S21.imag
#S21tofit = np.concatenate((it,qt))
#
#p0t = (f0t,Qrt,Qct,0,0)
#boundst = ([resdict['fine']['freqs'].min().value,1e3,1e3,-np.inf,-1],[resdict['fine']['freqs'].max().value,1e6,1e8,np.inf,1])
#res_popt,res_pcov = curve_fit(fit_S21_nonlinear,freqs.value,S21tofit,p0=p0t,bounds=boundst)
#plt.plot(freqs,np.abs(S21_nonlinear(freqs.value,*res_popt)),'c-',label='nonlinear S21 fit')
#plt.xlabel('freq (Hz)')
#plt.ylabel('|S21|')
#plt.legend(loc='lower left')
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
##
##    
