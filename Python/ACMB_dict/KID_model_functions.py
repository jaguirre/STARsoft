# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 22:18:24 2018

@author: Alyssa
"""

import numpy as np
from astropy import units as u
from astropy import constants as c
from scipy.special import kv,iv
from astropy.modeling.blackbody import blackbody_lambda,blackbody_nu
from scipy.integrate import trapz


N0_Al = 1.72e10*np.power(u.micron,-3)*np.power(u.eV,-1)

#%%

MB_const = 1.76
def delta0(Tc):
    if not hasattr(Tc,'unit'): Tc*=u.K

    delta = (MB_const*c.k_B*Tc).to(u.J)

    return delta

#%%
def n_gamma(f,T_BB):
    if not hasattr(T_BB,'unit'): T_BB*=u.K
    T_BB = T_BB.to(u.K)
    
    ex = (c.h*f/(c.k_B*T_BB)).to(u.dimensionless_unscaled)
    n_gamma = np.power(np.exp(ex)-1,-1)
    return n_gamma

#%%
filt_k,filt_t = np.loadtxt('BandpassFilter.txt',skiprows=1,delimiter=',',unpack=True)
indx = np.argsort(filt_k)
freq_filt = (filt_k[indx]*u.k).to(u.THz,equivalencies=u.spectral())
trans_filt = filt_t[indx]
from scipy.interpolate import CubicSpline
filt = CubicSpline(freq_filt.value,trans_filt) # function that takes in a freq in THz, returns a transmission value at that freq

def TBB_to_Pinc(TBB,trans=1):
    f = np.linspace(freq_filt.min(),freq_filt.max(),1000) # integrate over the range where we have filter transmission data
    B_nu = blackbody_nu(f,TBB).to(u.W/u.m**2/u.steradian/u.Hz)
    P_nu = B_nu*np.power(f.to(u.m,equivalencies=u.spectral()),2)*u.steradian
    tr = trans*filt(f.to(u.THz))*u.dimensionless_unscaled
    
    P = trapz(P_nu*tr,x=f).to(u.pW)
    
    return P

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
    S_1 = (2/np.pi)*np.sqrt((2*delta*beta)/(np.pi))*np.sinh(arg)*kv(0,arg)
    S_1 = S_1.to(u.dimensionless_unscaled)
    
    return S_1

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

    S_2 = 1 + np.sqrt((2*delta*beta)/(np.pi))*np.exp(-arg)*iv(0,arg)
    S_2 = S_2.to(u.dimensionless_unscaled)
    
    return S_2

#%%
def nth(Tstage,Tc,N0=N0_Al):
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc = Tc*u.K
    Tc = Tc.to(u.K)
    
    delta = delta0(Tc)
    n_th = (2*N0*(np.sqrt(2*np.pi*c.k_B*Tstage*delta)).to(u.eV)*np.exp(-delta/(c.k_B*Tstage))).to(np.power(u.micron,-3))
    
    return n_th

#%%
def nqp(Tstage,Tc,TBB,V,n_star,tau_max,eta_pb,trans=1):
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc*=u.K
    Tc = Tc.to(u.K)
    
    if not hasattr(TBB,'unit'): TBB*=u.K
    TBB = TBB.to(u.K)
    
    if not hasattr(V,'unit'): V = V*np.power(u.micron,3)
    V = V.to(np.power(u.micron,3))
    
    if not hasattr(tau_max,'unit'): tau_max*=u.microsecond
    tau_max = tau_max.to(u.microsecond)
    
    if not hasattr(n_star,'unit'): n_star*=np.power(u.micron,-3)
    n_star = n_star.to(np.power(u.micron,-3))
    

    eta_pb*=u.dimensionless_unscaled
    
    trans*=u.dimensionless_unscaled
    
    delta = delta0(Tc)
    
    Pinc = TBB_to_Pinc(TBB,trans)
    
    n_th = nth(Tstage,Tc)
    
    n_qp = (-n_star + ((n_star+n_th)**2 + (2*n_star*eta_pb*Pinc*tau_max)/(delta*V))**0.5).to(np.power(u.micron,-3))
    
    return n_qp

#%%
def tauqp(n_qp,tau_max,n_star):
    if not hasattr(n_qp,'unit'): n_qp*=np.power(u.micron,-3)
    n_qp = n_qp.to(np.power(u.micron,-3))
    
    if not hasattr(tau_max,'unit'): tau_max*=u.microsecond
    tau_max = tau_max.to(u.microsecond)
    
    if not hasattr(n_star,'unit'): n_star*=np.power(u.micron,-3)
    n_star = n_star.to(np.power(u.micron,-3))
    
    tau_qp = (tau_max/(1+n_qp/n_star)).to(u.microsecond)
    
    return tau_qp

#%%
def gammath(tau_max,n_star,V,Tstage,Tc):
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc*=u.K
    Tc = Tc.to(u.K)
    
    if not hasattr(V,'unit'): V = V*np.power(u.micron,3)
    V = V.to(np.power(u.micron,3))


    n_th = nth(Tstage,Tc)
    
    gamma_th = (((n_th*V)/2)*(np.power(tau_max,-1)+np.power(tauqp(n_th,tau_max,n_star),-1))).to(np.power(u.microsecond,-1))
    
    return gamma_th

#%%
def gammar(Tstage,Tc,TBB,V,n_star,tau_max,eta_pb,trans=1):
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc*=u.K
    Tc = Tc.to(u.K)
    
    if not hasattr(TBB,'unit'): TBB*=u.K
    TBB = TBB.to(u.K)
    
    if not hasattr(V,'unit'): V = V*np.power(u.micron,3)
    V = V.to(np.power(u.micron,3))

    if not hasattr(tau_max,'unit'): tau_max*=np.power(u.micron,-3)
    tau_max = tau_max.to(u.microsecond)
    
    if not hasattr(n_star,'unit'): n_star*=np.power(u.micron,-3)
    n_star = n_star.to(np.power(u.micron,-3))

    Pinc = TBB_to_Pinc(TBB,trans)
    
    n_qp = nqp(Tstage,Tc,TBB,V,n_star,tau_max,eta_pb,trans)
    tau_qp = tauqp(n_qp,tau_max,n_star)
    
    gamma_r = (((n_qp*V)/2)*(np.power(tau_max,-1)+np.power(tau_qp,-1))).to(np.power(u.microsecond,-1))
    return gamma_r

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
def Q_inv_MB(alpha,Tc,Tstage,f,nqp,N0=N0_Al):
    delta = delta0(Tc)
    
    Q_inv_MB = nqp*(alpha*S1(Tstage,f,Tc))/(2*N0*delta)
    Q_inv_MB = Q_inv_MB.to(u.dimensionless_unscaled)
    
    return Q_inv_MB

#%%
def Q_inv_MB_fit(data,alpha,Tc,Q_0):
    Tstage = data[0]
    f = data[1]
    
    nqp = nth(Tstage,Tc)
    Q_inv_MB_fit = Q_inv_MB(alpha,Tc,Tstage,f,nqp)
    
    return Q_inv_MB_fit

#%%
def x_opt(alpha,Tstage,f,Tc,tau_max,n_star,Pinc,eta_pb,N0=N0_Al):
    S_2 = S2(Tstage,f,Tc)
    delta = delta0(Tc)
    P_abs = eta_opt*Pinc
    
    x_opt = (-alpha*S_2/(4*N0*delta)*(n_star*(1.+(2*eta_pb*P_abs*tau_max)/(delta*V*n_star))**0.5 - n_star)).to(u.dimensionless_unscaled)
    
    return x_opt

#%%
def Sxx_tot(alpha,Tstage,f,Tc,tau_max,n_star,Pinc,V,eta_pb,nu_opt,eta_opt,N0=N0_Al,n_gamma=0):
    S_2 = S2(Tstage,f,Tc)
    P_abs = eta_opt*Pinc
    n_qp = nqp(Tstage,Tc,P_abs,V,n_star,tau_max,eta_pb)
    tau_qp = tauqp(n_qp,tau_max,n_star)
    gamma_th = gammath(tau_max,n_star,V,Tstage,Tc)
    gamma_r = gammar(Tstage,Tc,P_abs,V,n_star,tau_max,eta_pb)
    delta = delta0(Tc)
    N0 = N0.to(np.power(u.J,-1)*np.power(u.micron,-3))
    #n_gamma = 1/e^(h*f/kT -1)
       
    Sxx_tot = ((alpha*S_2/(4*N0*delta))**2 * ((((eta_pb*tau_qp/(delta*V))**2)*2*c.h*f*P_abs*(1+n_gamma))+((4*tau_qp**2/V**2)*(gamma_th+gamma_r)))).to(np.power(u.Hz,-1))
    return Sxx_tot

#%%
#%%

#test case:
alpha = 0.8*u.dimensionless_unscaled
Tstage = 0.215*u.K
f = 300*u.MHz
Tc = 1.4*u.K
tau_max = 35*u.microsecond
n_star = 1318*(np.power(u.micron,-3))
Pinc = 0.5*u.pW
V = 76*np.power(u.micron,3)
eta_pb = 0.57
nu_opt = (350*u.micron).to(u.GHz,equivalencies=u.spectral())
eta_opt = 0.17*u.dimensionless_unscaled
N0=N0_Al
n_gamma=0

Pinc = np.linspace(0.01,2.5,20)*u.pW
sxxtest = Sxx_tot(alpha,Tstage,f,Tc,tau_max,n_star,Pinc,V,eta_pb,nu_opt,eta_opt,N0,n_gamma)
xtest = x_opt(alpha,Tstage,f,Tc,tau_max,n_star,Pinc,eta_pb)
#%%
# individual device fitting

# import dark data for individual resonator
cd011_res3 = np.transpose(np.loadtxt('cd011_res3.csv',delimiter=',',skiprows=1))
T_stage,xavg,xerr,Qravg,Qrerr,Sxx_avg,Sxx_err = np.loadtxt('cd011_res3.csv',delimiter=',',skiprows=1,unpack=True)
T_stage*=u.K
xavg*=u.dimensionless_unscaled
xerr*=u.dimensionless_unscaled
Qravg*=u.dimensionless_unscaled
Qrerr*=u.dimensionless_unscaled
Sxx_avg*=np.power(u.Hz,-1)
Sxx_err*=np.power(u.Hz,-1)
#for el in [xavg,xerr,Qravg,Qrerr]: el[0]*=u.dimensionless_unscaled
#for el in [Sxx_avg,Sxx_err]: el*=(np.power(u.Hz,-1))



# fit dark data: simultaneous fit of x and Qinv vs T_stage --> free parameters alpha, Tc, dx0, Qinv0


# import optical data for individual resonator
T_BB,xavg,xerr,Qravg,Qrerr,Sxx_avg,Sxx_err = np.loadtxt('cd010_res1.csv',delimiter=',',skiprows=1,unpack=True)
for el in [xavg,xerr,Qravg,Qrerr]: el*=u.dimensionless_unscaled
for el in [Sxx_avg,Sxx_err]: el*=(np.power(u.Hz,-1))






