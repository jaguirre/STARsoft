# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 13:59:17 2018

@author: Alyssa
"""

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

#Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt,trans=1,N0=N0_Al):
T_BB = np.linspace(.1,25,num=100)*u.K

alpha=0.8
f=330*u.MHz
Tstage = .215*u.K
Tc=1.35
V=76*(38/40)*.8*np.power(u.micron,3)
eta_pb = 0.57
nu_opt = (350*u.micron).to(u.GHz,equivalencies=u.spectral())
eta_opt = 1

n_star = 100*np.power(u.micron,-3)
tau_max = 40*u.microsecond
plt.plot(TBB_to_Pinc(T_BB),Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt),'r-')

tau_max = 100*u.microsecond
plt.plot(TBB_to_Pinc(T_BB),Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt),'b-')

tau_max = 1000*u.microsecond
plt.plot(TBB_to_Pinc(T_BB),Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt),'g-')

tau_max = 75*u.microsecond
plt.plot(TBB_to_Pinc(T_BB),Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt),'m-')

n_star = 1000*np.power(u.micron,-3)
tau_max = 40*u.microsecond
plt.plot(TBB_to_Pinc(T_BB),Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt),'r:')

eta_opt = 0.1
plt.plot(TBB_to_Pinc(T_BB),Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt),'r--')
