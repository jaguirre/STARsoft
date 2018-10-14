# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 11:28:33 2018

@author: Alyssa
"""

import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c
from dictionary_functions import *
import KID_model_functions as kids

#%%
from astropy import units as u

def TBB_to_Pinc(TBB,trans=1):
    '''function to convert blackbody temperature to incident optical power for the cryogenic blackbody
        in the Red cryostat at Caltech. If a transmission mask is in place, give its transmission efficiency 
        as a decimal value using the trans argument (1=no mask, 0=totally dark) '''
        
    bbtemps,bbpowers = np.loadtxt('Optical_power_vs_T_BB.csv',delimiter=',',unpack=True)
    bbtemps*=u.K
    bbpowers*=u.pW
    
    if not hasattr(TBB,'unit'): TBB*=u.K #assume TBB is in K if no units given
    else: TBB = TBB.to(u.K)
    
    if TBB.isscalar: Pinc = bbpowers[max(np.where(bbtemps < TBB)[0]+1)] # scalars are special
    else: Pinc = [bbpowers[max(np.where(bbtemps < temp)[0]+1)] for temp in TBB]
    
    Pinc*=trans
    return Pinc

#%%

''' hdf5 files packed up on hoverboard with dd.io.save '''
cd010 = dd.io.load('CD010_data.hdf5')
cd011 = dd.io.load('CD011_data.hdf5')
cd012 = dd.io.load('CD012_data.hdf5')
cd013 = dd.io.load('cd013_data.hdf5')

''' Optical transmssion for each cooldown '''
cd010['T_opt'] = 1 # fully open, no mask
cd011['T_opt'] = 0 # fully dark
cd012['T_opt'] = 0.03 # open with 3% mask in place
cd013['T_opt'] = 1 # fully open, no mask


cdlist = [cd010,cd011,cd012,cd013]

#%%

''' Resonance frequencies of resonators we studied for this device: '''
res0 = 324.56*u.MHz
res1 = 331.6*u.MHz
res2 = 336.84*u.MHz
res3 = 342.22*u.MHz
res4 = 352.03*u.MHz
res5 = 355.61*u.MHz
res6 = 361.85*u.MHz
res7 = 363.59*u.MHz
resonators=[res0,res1,res2,res3,res4,res5,res6,res7]


''' We'll use the range df=f0/Qrange to determine which resonator is which
    [not the most elegant solution if there are lots of tightly-packed resonators and collisions] '''
Qrange = 300
df = [f0/Qrange for f0 in resonators]



for cool in cdlist:
    
    ''' Compare each resonance frequency to the known list of resonators, 
        use the value to index each scan by resonator number '''
    for freq in list(gen_dict_extract('f0',cool)):
        path = dictwhere(cool,'f0',freq)
        for ind,res in enumerate(resonators):
            if abs(res-freq) < df[ind]:
                for sweep,scan in path: # hard-coding in that each freq will be indexed by cd[sweep][scan]
                    cool[sweep][scan]['res'] = ind    
    
    ''' Find the maximum resonance frequency for each resonator within a cooldown,
        set the maximum to be f00 for that resonator,cooldown '''
    f00 = np.zeros(len(resonators))*u.MHz
    for n,r in enumerate(resonators):
        try: 
            f0n = max(dictget(cool,dictwhere(cool,'res',n),'f0'))
            f00[n] = f0n
            
        # if a resonator doesn't appear in this cooldown, set its f00 = -1
        except ValueError: 
            f00[n] = -1*u.MHz
    cool['f00'] = f00
    
    ''' Use the f00 values for each resonator,cooldown to calculate x = (f-f00)/f
        for each scan in the cooldown '''
    for freq in list(gen_dict_extract('f0',cool)):
        path2 = dictwhere(cool,'f0',freq)
        reso = dictget(cool,path2,'res')
        x = ((freq-f00[reso])/f00[reso]).to(u.dimensionless_unscaled)
        for sweep2,scan2 in path2: # hard-coding in that each freq will be indexed by cd[sweep][scan]
            cool[sweep2][scan2]['x'] = x
            
    del path,path2        
#%%            
f0 = plt.figure(100)
p0 = f0.add_subplot(111)
   
#%%

colors = ['r','orange','coral','y','g','limegreen','b','purple']

''' Plot: x and Sxx vs P_inc for an individual resonator for cd012, for the lowest value of LB_atten
'''
d = cd012

#plt.close('all')
save = False

for reso in np.arange(6,7):#len(resonators)):
    
    conds = [['res',reso],['T_stage',215*u.mK]]
    resodict = redict(d,dictcut(d,conds)) # dictionary containing only the elements for this resonator    
    
    if resodict!={}:
        print(reso)
        pows = list(gen_dict_extract('P_BB',resodict)) 
        upows = list(set(pows)) 
        upows.sort()
        upows.remove(upows[2]) ### CD012, something weird happened at this point (mis-recorded LB power?)
        
        f1 = plt.figure(reso)
        p1 = f1.add_subplot(111)
        
        f2 = plt.figure('-'+str(reso))
        p2 = f2.add_subplot(111)
    
        pavg = []
        xavg = []
        x_err = []
        sxxavg = []
        sxxerr = []
        # take the min value of LB_atten for each resonator,T_stage,optical power config
        for popt in upows: 
            paths = dictwhere(resodict,'P_BB',popt)
            tdict = redict(resodict,paths)
            tonepows = list(gen_dict_extract('LB_atten',tdict))
            db3 = min(tonepows) # hard-coding that the lowest attenuation we'll use is always ~3 dB below bifurcation (ish)
        
            paths2 = dictwhere(tdict,'LB_atten',db3)
            tdict2 = redict(tdict,paths2)
            
            xlist = list(gen_dict_extract('x',tdict2))
            xlist = [x.to(u.dimensionless_unscaled).value for x in xlist] # this is dumb but necessary to plot
            xavg.append(np.mean(xlist))
            x_err.append(np.std(xlist)/np.sqrt(len(xlist)))

            sxxlist = list(gen_dict_extract('CSD_bin_subtracted_avg',tdict2))
            sxxlist = [sxx.to(np.power(u.Hz,-1)).value for sxx in sxxlist] # this is dumb but necessary to plot
            sxxavg.append(np.mean(sxxlist))
            sxxerr.append(np.std(sxxlist)/np.sqrt(len(sxxlist)))
            
            plist = list(gen_dict_extract('P_BB',tdict2)) # should be all equal to popt!
            plist = [p.to(u.pW).value for p in plist] # this is dumb but necessary to plot
            pavg.append(plist[0])
        
            p1.plot(plist,xlist,'k.',alpha=0.75) # plot the individual points
            p2.plot(plist,sxxlist,'k.',alpha=0.75)
            
        p1.errorbar(pavg,xavg,yerr=x_err,color='b',fmt='o',capsize=2.0,ecolor='b')
        p2.errorbar(pavg,sxxavg,yerr=sxxerr,color='b',fmt='o',capsize=2.0,ecolor='b')
        
        resp = np.polyfit(pavg,xavg,1)
        print(resp[0])
        
        p1.plot(pavg,np.polyval(resp,pavg),'g-')
        color=colors[reso]
        p0.errorbar(pavg,xavg,yerr=x_err,color=color,fmt='s',capsize=2.0,ecolor=color)
        p0.errorbar(list(d['T_opt']*np.asarray(pavg)),xavg,yerr=x_err,color=color,fmt='*',capsize=2.0,ecolor=color)

        
        p1.grid(color='gray',linestyle='--',alpha=0.75)
        p1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        p1.tick_params(direction='in')
        p1.set_title('CD012, res ' + str(reso),fontsize=14)
        p1.set_xlabel(r'P$_{BB}$ (pW)',fontsize=14)
        p1.set_xlim([0.05,0.55])
        p1.set_ylabel(r'$\delta$f/f',fontsize=14)
        p1.set_ylim([-6e-6,0])
        p1.tick_params(axis='both', which='major', labelsize=12)
        p1.text(0.1,-5e-6, 'R = '+ "{:.2e}".format(resp[0]) + ' /pW',fontsize=12)
        if save==True: f1.savefig('results/cd012_res'+str(reso)+'Pinc_vs_x.png')
        
        p2.grid(color='gray',linestyle='--',alpha=0.75)
        p2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        p2.tick_params(direction='in')
        p2.set_title('CD012, res ' + str(reso),fontsize=14)
        p2.set_xlabel(r'P$_{BB}$ (pW)',fontsize=14)
        p2.set_xlim([0.05,0.55])
        p2.set_ylabel(r'$S_{xx}$ white noise (1/Hz)',fontsize=14)
#        p2.set_ylim([-6e-6,0])
        p2.tick_params(axis='both', which='major', labelsize=12)
        if save==True: f2.savefig('results/cd012_res'+str(reso)+'Pinc_vs_sxx.png')

#%%
''' Plot: x and Sxx vs P_inc for an individual resonator for cd010, for the lowest value of LB_atten
'''

d = cd010

#plt.close('all')
save = False

for reso in np.arange(0,len(resonators)):
    
    conds = [['res',reso],['T_stage',215*u.mK]]
    resodict = redict(d,dictcut(d,conds)) # dictionary containing only the elements for this resonator    
    
    if resodict!={}:
        pows = list(gen_dict_extract('P_inc',resodict)) 
        upows = list(set(pows)) 
        upows.sort()
        #upows.remove(upows[2]) ### CD012, something weird happened at this point (mis-recorded LB power?)
        
        f1 = plt.figure(reso)
        p1 = f1.add_subplot(111)
        
        f2 = plt.figure('-'+str(reso))
        p2 = f2.add_subplot(111)
    
        pavg = []
        xavg = []
        x_err = []
        sxxavg = []
        sxxerr = []
        # take the min value of LB_atten for each resonator,T_stage,optical power config
        for popt in upows: 
            paths = dictwhere(resodict,'P_inc',popt)
            tdict = redict(resodict,paths)
            tonepows = list(gen_dict_extract('LB_atten',tdict))
            db3 = min(tonepows) # hard-coding that the lowest attenuation we'll use is always ~3 dB below bifurcation (ish)
        
            paths2 = dictwhere(tdict,'LB_atten',db3)
            tdict2 = redict(tdict,paths2)
            
            xlist = list(gen_dict_extract('x',tdict2))
            xlist = [x.to(u.dimensionless_unscaled).value for x in xlist] # this is dumb but necessary to plot
            xavg.append(np.mean(xlist))
            x_err.append(np.std(xlist)/np.sqrt(len(xlist)))

            sxxlist = list(gen_dict_extract('CSD_bin_subtracted_avg',tdict2))
            sxxlist = [sxx.to(np.power(u.Hz,-1)).value for sxx in sxxlist] # this is dumb but necessary to plot
            sxxavg.append(np.mean(sxxlist))
            sxxerr.append(np.std(sxxlist)/np.sqrt(len(sxxlist)))
            
            plist = list(gen_dict_extract('P_inc',tdict2)) # should be all equal to popt!
            plist = [p.to(u.pW).value for p in plist] # this is dumb but necessary to plot
            pavg.append(plist[0])
        
            p1.plot(plist,xlist,'k.',alpha=0.75) # plot the individual points
            p2.plot(plist,sxxlist,'k.',alpha=0.75)
            
        p1.errorbar(pavg,xavg,yerr=x_err,color='b',fmt='o',capsize=2.0,ecolor='b')
        p2.errorbar(pavg,sxxavg,yerr=sxxerr,color='b',fmt='o',capsize=2.0,ecolor='b')
        
        resp = np.polyfit(pavg,xavg,1)
        print(resp[0])
        
        p1.plot(pavg,np.polyval(resp,pavg),'g-')
        
        color=colors[reso]
        p0.errorbar(pavg,xavg,yerr=x_err,color=color,fmt='o',capsize=2.0,ecolor=color)

        
        p1.grid(color='gray',linestyle='--',alpha=0.75)
        p1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        p1.tick_params(direction='in')
        p1.set_title('CD010, res ' + str(reso),fontsize=14)
        p1.set_xlabel(r'P$_{inc}$ (pW)',fontsize=14)
        p1.set_xlim([0.1,2.5])
        p1.set_ylabel(r'$\delta$f/f',fontsize=14)
        p1.set_ylim([-4.25e-4,1e-5])
        p1.tick_params(axis='both', which='major', labelsize=12)
        p1.text(0.25,-4e-4, 'R = '+ "{:.2e}".format(resp[0]) + ' /pW',fontsize=12)
        if save==True: f1.savefig('results/cd010_res'+str(reso)+'Pinc_vs_x.png')
        
        p2.grid(color='gray',linestyle='--',alpha=0.75)
        p2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        p2.tick_params(direction='in')
        p2.set_title('CD010, res ' + str(reso),fontsize=14)
        p2.set_xlabel(r'P$_{inc}$ (pW)',fontsize=14)
        p2.set_xlim([0.1,2.5])
        p2.set_ylabel(r'$S_{xx}$ white noise (1/Hz)',fontsize=14)
        p2.set_ylim([0.4e-16,2.2e-16])
        p2.tick_params(axis='both', which='major', labelsize=12)
        if save==True: f2.savefig('results/cd010_res'+str(reso)+'Pinc_vs_sxx.png')
#%%
d = cd013

#plt.close('all')
save = False

for reso in np.arange(0,len(resonators)):
    
    conds = [['res',reso],['T_stage',215*u.mK]]
    resodict = redict(d,dictcut(d,conds)) # dictionary containing only the elements for this resonator    
    
    if resodict!={}:
        pows = list(gen_dict_extract('P_inc',resodict)) 
        upows = list(set(pows)) 
        upows.sort()
        #upows.remove(upows[2]) ### CD012, something weird happened at this point (mis-recorded LB power?)
        
        f1 = plt.figure(reso)
        p1 = f1.add_subplot(111)
        
        f2 = plt.figure('-'+str(reso))
        p2 = f2.add_subplot(111)
    
        pavg = []
        xavg = []
        x_err = []
        sxxavg = []
        sxxerr = []
        # take the min value of LB_atten for each resonator,T_stage,optical power config
        for popt in upows: 
            paths = dictwhere(resodict,'P_inc',popt)
            tdict = redict(resodict,paths)
            tonepows = list(gen_dict_extract('LB_atten',tdict))
            db3 = min(tonepows) # hard-coding that the lowest attenuation we'll use is always ~3 dB below bifurcation (ish)
        
            paths2 = dictwhere(tdict,'LB_atten',db3)
            tdict2 = redict(tdict,paths2)
            
            xlist = list(gen_dict_extract('x',tdict2))
            xlist = [x.to(u.dimensionless_unscaled).value for x in xlist] # this is dumb but necessary to plot
            xavg.append(np.mean(xlist))
            x_err.append(np.std(xlist)/np.sqrt(len(xlist)))

            sxxlist = list(gen_dict_extract('CSD_bin_subtracted_avg',tdict2))
            sxxlist = [sxx.to(np.power(u.Hz,-1)).value for sxx in sxxlist] # this is dumb but necessary to plot
            sxxavg.append(np.mean(sxxlist))
            sxxerr.append(np.std(sxxlist)/np.sqrt(len(sxxlist)))
            
            plist = list(gen_dict_extract('P_inc',tdict2)) # should be all equal to popt!
            plist = [p.to(u.pW).value for p in plist] # this is dumb but necessary to plot
            pavg.append(plist[0])
        
            p1.plot(plist,xlist,'k.',alpha=0.75) # plot the individual points
            p2.plot(plist,sxxlist,'k.',alpha=0.75)
            
        p1.errorbar(pavg,xavg,yerr=x_err,color='b',fmt='o',capsize=2.0,ecolor='b')
        p2.errorbar(pavg,sxxavg,yerr=sxxerr,color='b',fmt='o',capsize=2.0,ecolor='b')
        
        resp = np.polyfit(pavg,xavg,1)
        print(resp[0])
        
        p1.plot(pavg,np.polyval(resp,pavg),'g-')
        
        color=colors[reso]
        p0.errorbar(pavg,xavg,yerr=x_err,color=color,fmt='^',capsize=2.0,ecolor=color)

        
        p1.grid(color='gray',linestyle='--',alpha=0.75)
        p1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        p1.tick_params(direction='in')
        p1.set_title('CD010, res ' + str(reso),fontsize=14)
        p1.set_xlabel(r'P$_{inc}$ (pW)',fontsize=14)
        p1.set_xlim([0.1,2.5])
        p1.set_ylabel(r'$\delta$f/f',fontsize=14)
        p1.set_ylim([-4.25e-4,1e-5])
        p1.tick_params(axis='both', which='major', labelsize=12)
        p1.text(0.25,-4e-4, 'R = '+ "{:.2e}".format(resp[0]) + ' /pW',fontsize=12)
        if save==True: f1.savefig('results/cd010_res'+str(reso)+'Pinc_vs_x.png')
        
        p2.grid(color='gray',linestyle='--',alpha=0.75)
        p2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        p2.tick_params(direction='in')
        p2.set_title('CD010, res ' + str(reso),fontsize=14)
        p2.set_xlabel(r'P$_{inc}$ (pW)',fontsize=14)
        p2.set_xlim([0.1,2.5])
        p2.set_ylabel(r'$S_{xx}$ white noise (1/Hz)',fontsize=14)
        p2.set_ylim([0.4e-16,2.2e-16])
        p2.tick_params(axis='both', which='major', labelsize=12)
        if save==True: f2.savefig('results/cd010_res'+str(reso)+'Pinc_vs_sxx.png')

        
        