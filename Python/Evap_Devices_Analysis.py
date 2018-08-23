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
        reso = dictget(cool,path,'res')
        x = ((freq-f00[reso])/f00[reso]).to(u.dimensionless_unscaled)
        for sweep2,scan2 in path2: # hard-coding in that each freq will be indexed by cd[sweep][scan]
            cool[sweep2][scan2]['x'] = x
            
    ''' For each scan, translate T_BB into P_inc using the calculated table 
        and optical transmission for that cooldown '''
    trans = cool['T_opt']
    for temp in list(gen_dict_extract('T_BB',cool)):
        p_inc = TBB_to_Pinc(temp,trans)
        path3 = dictwhere(cool,'T_BB',temp)
        for sweep3,scan3 in path3: # hard-coding in that each bb temp will be indexed by cd[sweep][scan]
            cool[sweep3][scan3]['P_inc'] = p_inc
            
f1 = plt.figure(1)
p1 = f1.add_subplot(111)

pows = list(gen_dict_extract('P_inc',cd011))
for ind,power in enumerate(pows):
    x = list(gen_dict_extract('x',cd011))[ind]
    p1.plot(power,x,'k.')            
            
