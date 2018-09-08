#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 18:46:25 2017

@author: jaguirre
"""

from astropy import units as u
from astropy import constants as c
#from astropy.analytic_functions import blackbody_lambda
#from astropy.analytic_functions import blackbody_nu
from astropy.modeling.blackbody import blackbody_lambda,blackbody_nu
from astropy.modeling import models
from scipy.integrate import trapz

import numpy as np
import matplotlib.pylab as plt

def Qprint(quantity,sigfig=3,unitstyle='s'):
    """ wrap up the syntax for printing astropy Quantity objects in a pretty 
    way. Other options for unitstyle are latex and ... ?  """
    sf = str(sigfig) # Need to make a zero padded thing to cope with the case
    # of many sigfigs
    fmtstr = '{0.value:0.00'+str(sigfig)+'g} {0.unit:'+unitstyle+'}'
    x = fmtstr.format(quantity) 
    return x

def B_lambda(wavelength,T):
    B = blackbody_lambda(wavelength,T)
    B = B.to(u.W/u.m**3/u.steradian)
    return B

def B_nu(frequency,T):
    B = blackbody_nu(frequency,T)
    B = B.to(u.W/u.m**2/u.steradian/u.Hz)
    return B

def xBB(nu,T):
    return (c.h*nu/(c.k_B*T)).to(u.dimensionless_unscaled)

def n_nu(nu,T):
    return np.power(np.exp(xBB(nu,T)) - 1, -1)

#def B_nu(nu,T):
#    n = n_nu(nu,T) #np.power(np.exp(c.h*nu/(c.k_B * T)) - 1, -1)
#    B = 2.*c.h*np.power(nu,3)/np.power(c.c,2) * n / u.sr
#    B = B.to(u.W/u.m**2/u.steradian/u.Hz)
#    return B

#%%
def deriv(f,x):
    df = np.gradient(f)
    # allows for the possibility that dx is not constant.  who knows how accurate?
    dx = np.gradient(x)
    return df/dx


#%%
def dB_nu_dT(nu,T):
    dB = (2.*c.h*np.power(nu,3)/np.power(c.c,2) * np.power(n_nu(nu,T),2) * np.exp(xBB(nu,T)) * 
    xBB(nu,T)/T).to(u.W/u.m**2/u.Hz/u.K)/u.sr
    return dB
#%%

Ts = np.arange(4,21)*u.K#,100,500,1000,2000])*u.K
nu = np.linspace(1e-3,10,num=1e6)*u.THz
lmbda = (c.c/nu).to(u.m)

nu0 = (850*u.GHz).to(u.THz)
fwhm = 0.12*nu0
sigma = fwhm/2.355
g = models.Gaussian1D(amplitude=0.9, mean=nu0.value, stddev=sigma.value)
filt = g(nu.value)
dnu = trapz(filt,x=nu)

P = np.zeros(len(Ts))*u.W

plt.figure(1)
plt.clf()
for i,T in enumerate(Ts):
    B = B_nu(nu,T)
    P_nu = B*np.power(lmbda,2)*u.steradian
    # only for plotting
    norm = (B_nu(nu,Ts.max())*np.power(lmbda,2)*u.steradian).max()
    plt.semilogy(nu,P_nu/norm,label='T = '+str(T))
    P[i] = trapz(P_nu*filt,x=nu)
    print(T,'  ',Qprint(P[i].to(u.pW)))
plt.plot(nu,filt,'k')
plt.xlim([0,2])
plt.xlabel('Frequency [THz]')
plt.ylim([1e-10,1])
plt.ylabel(r'$\propto B_\nu(T,\nu) c^2/\nu^2$')
plt.legend(loc='lower left')


#%%
plt.figure(2)
plt.clf()
plt.semilogy(Ts,P.to(u.pW),label=r'$\int B_\nu(T,\nu) \lambda^2 W(\nu) d\nu$' )
plt.semilogy(Ts,(2.*c.k_B*Ts*dnu).to(u.pW),label=r'$2 k T \int W(\nu) d\nu $')
plt.legend(loc='lower right')
plt.ylabel('Power [pW]')
plt.xlabel('Temperature [K]')
plt.grid()
plt.savefig('PvsT.png')
#%%
plt.figure(3)
plt.clf()
dlnPdT = np.zeros(len(Ts))/u.mK
for i,T in enumerate(Ts):
    dlnPdT[i] = trapz(dB_nu_dT(nu,T)*np.power(lmbda,2)*u.sr*filt,x=nu)/trapz(B_nu(nu,T)*np.power(lmbda,2)*u.sr*filt,x=nu)
    #plt.plot(nu,dlnPdT,label='T = '+str(T))

dlnPdT.to(1/u.mK)

#%%
plt.plot(Ts,dlnPdT,label=r'Planck at 350$\mu$m')
plt.plot(Ts,(1/Ts).to(1/u.mK),label=r'RJ at 350$\mu$m')
plt.ylabel(r'$\frac{1}{P}\frac{dP}{dT}$ [mK$^{-1}$]')
plt.xlabel('Temperature [K]')
plt.legend(loc='upper right')
plt.savefig('dlnPdTvsT.png')