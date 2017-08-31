#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 16:37:07 2017

@author: jaguirre
"""

import numpy as np
from astropy.cosmology import Planck15 as cosmo
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c


#%%

def phi_bethermin11(L_IR,z):
    phi_star = 1.
    L_star = 1.e11
    phi = phi_star * np.power(L_IR/L_star,-1)
    return phi

def phi_somebody_else(L_IR,z):
    phi_star = 1.
    L_star = 1.e11
    phi = phi_star * np.power(L_IR/L_star,-1) * np.exp(-L_IR/1e11)
    return phi

def Sofz(luminosity_function,z,species='CII'):
    
    if species == 'CII':
        f = 0.1
    else:
        f = 1.
    L_IR = np.linspace(1e8,1e13)
    phi = f*luminosity_function(L_IR,z)
    S = np.trapz(phi, x=L_IR)
    return S

#%%
def dchidnu(z,lambda_rest):
    return lambda_rest * np.power(1+z,2) / cosmo.H(z)
    

#%%

L_IR = np.linspace(1e8,1e13)
plt.loglog(L_IR,phi_bethermin11(L_IR,0))
plt.loglog(L_IR,phi_somebody_else(L_IR,0))
print Sofz(phi_bethermin11,0.)
print Sofz(phi_somebody_else,0.)

#%%
z = np.linspace(0.5,3,num=100)
d_A = cosmo.angular_diameter_distance(z)
d_L = cosmo.luminosity_distance(z)
y = dchidnu(z,158*u.micron).to(u.Mpc/u.Hz)
#%%
plt.figure(1)
plt.clf()
#plt.plot(z,d_A)
#plt.plot(z,d_L)
plt.plot(z,y) 