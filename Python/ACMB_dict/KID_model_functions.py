# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 22:18:24 2018

@author: Alyssa
"""

import numpy as np
from astropy import units as u
from astropy import constants as c
from scipy.special import kv,iv

N0_Al = 1.72e10*np.power(u.micron,-3)*np.power(u.eV,-1)

#%%

MB_const = 1.76
def delta0(Tc):
    if not hasattr(Tc,'unit'): Tc*=u.K
    delta = MB_const*c.k_B*Tc
    return delta

#%%
def S1(Tstage,f,Tc):
    if not hasattr(Tstage,'unit'):
        if Tstage > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
    
    if not hasattr(Tc,'unit'): Tc = Tc*u.K
    Tc = Tc.to(u.K)

    if not hasattr(f,'unit'): 
        if f > 1: f = f*u.MHz
        else: f = f*u.GHz
    f = f.to(u.MHz)

    beta = 1/(c.k_B*Tstage)
    
    delta = delta0(Tc)
    
    arg = (0.5*beta*c.h*f).to(u.dimensionless_unscaled).value
    S1 = (2/np.pi)*np.sqrt((2*delta*beta)/(np.pi))*np.sinh(arg)*kv(0,arg)
    S1 = S1.to(u.dimensionless_unscaled)
    
    return S1
#%%
def S2(Tstage,f,Tc):
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
    
    if not hasattr(Tc,'unit'): Tc = Tc*u.K
    Tc = Tc.to(u.K)

    if not hasattr(f,'unit'): 
        if f > 1: f = f*u.MHz
        else: f = f*u.GHz
    f = f.to(u.MHz)

    beta = 1/(c.k_B*Tstage)
    
    delta = delta0(Tc)
    
    arg = (0.5*beta*c.h*f).to(u.dimensionless_unscaled).value

    S2 = 1 + np.sqrt((2*delta*beta)/(np.pi))*np.exp(-arg)*iv(0,arg)
    S2 = S2.to(u.dimensionless_unscaled)
    
    return S2
#%%
def nth(Tstage,Tc,N0=N0_Al):
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc = Tc*u.K
    Tc = Tc.to(u.K)
    
    delta = delta0(Tc)
    nth = (2*N0*(np.sqrt(2*np.pi*c.k_B*Tstage*delta)).to(u.eV)*np.exp(-delta/(c.k_B*Tstage))).to(np.power(u.micron,-3))
    
    return nth
#%%
def nqp(Tstage,Tc,Pabs,V,n_star,tau_max,eta_pb):
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc = Tc*u.K
    Tc = Tc.to(u.K)
    
    if not hasattr(Pabs,'unit'): Pabs = Pabs*u.pW
    Pabs = Pabs.to(u)
    
    if not hasattr(V,'unit'): V = V*np.power(u.micron,3)
    V = V.to()


#%%
def x_MB(alpha,Tc,Tstage,f,nqp,N0=N0_Al):
    delta = delta0(Tc)
    
    x_MB = -nqp*(alpha*S2(Tstage,f,Tc))/(4*N0*delta)
    x_MB = x_MB.to(u.dimensionless_unscaled)
    
    return x_MB
#%%
def x_MB_fit(data,alpha,Tc,dx):
    Tstage = data[0]
    f = data[1]
    
    nqp = nth(Tstage,Tc)

    x_MB_fit = x_MB(alpha,Tc,Tstage,f,nqp)+dx

    return x_MB_fit
#%%
def x_MB_fit_2(data,alpha,Tc):
    Tstage = data[0]
    f = data[1]
    
    nqp = nth(Tstage,Tc)

    x_MB_fit = x_MB(alpha,Tc,Tstage,f,nqp)

    return x_MB_fit



#%%
#test case:
Tstage = 0.2*u.K
Tc = 1.4*u.K
f = 200*u.MHz
nqp = nth(Tstage,Tc)
alpha = 0.8

    
    