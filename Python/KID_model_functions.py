# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 22:18:24 2018

@author: Alyssa
"""

import numpy as np
from astropy import units as u
from astropy import constants as c
import os
# This determines the actual location of the module file, so we can specify paths relative to it
path = os.path.abspath(__file__)
module_path = os.path.dirname(path)+os.sep

# The accepted literature value for the density of states in thin film Al. 
# It's set as default for all these functions so changing it here should change it everywhere.
N0_Al = 1.72e10*np.power(u.micron,-3)*np.power(u.eV,-1) 

# Electron-phonon interaction time in Al, from Kaplan++ 1976
tau_0_Al = 438e-9*u.second

#%%
''' Function to calculate the superconducting band gap delta0 as a function of critical temperature Tc 
    The parameter relating them, MB_const, could be adjusted from its typical Mattis-Bardeen value of 1.76 for non-Mattis-Bardeen superconductors
     Inputs -- Tc: Superconducting critical temperature in K
    Outputs -- delta: Superconducting bandgap energy delta in J '''
MB_const = 1.76
def delta0(T_c):
    # Get our units in order 
    if not hasattr(T_c,'unit'): T_c*=u.K
    T_c = T_c.to(u.K)
    # Calculate delta
#    print(T_c)
    delta = (MB_const*c.k_B*T_c).to(u.J)
#    delta = (MB_const*c.k_B*T_c)


    return delta

#%%
''' Function to calculate the photon occupation number n_gamma (also called n_0) of photons of frequency nu_opt,
    emitted from a blackbody source at temperature T_BB
     Inputs -- nu_opt: frequency of optical photons (in Hz or similar)
               T_BB: blackbody temperature in K (can be an array)
    Outputs -- n_gamma: photon occupation number (unitless)
'''
def ngamma(nu_opt,T_BB):
    # Get our units in order 
    if not hasattr(T_BB,'unit'): T_BB*=u.K
    T_BB = T_BB.to(u.K)
    
    # For convenience, calculate the exponent first
    ex = (c.h*nu_opt/(c.k_B*T_BB)).to(u.dimensionless_unscaled)
    
    # Calculate n_gamma
    n_gamma = np.power(np.exp(ex)-1,-1)
    return n_gamma

#%%
''' Function to calculate n_star, based on the electron-phonon interaction time and recombination constant
     Inputs -- tau_max: quasiparticle lifetime constant in microseconds
               Tc: Critical temperature of the superconductor in K
               tau_0: electron-phonon interaction time in seconds (default is reference value for Al)
               N0: Density of states for the superconductor material (default is thin-film Al value)
    Outputs -- quasiparticle number constant in microns^-3
'''
def nstar(tau_max,Tc,tau_0=tau_0_Al,N0=N0_Al):
    # Get our units in order
    if not hasattr(Tc,'unit'): Tc = Tc*u.K
    Tc = Tc.to(u.K)
    
    if not hasattr(tau_max,'unit') or tau_max.unit==u.dimensionless_unscaled: tau_max*=u.microsecond
    tau_max = tau_max.to(u.microsecond)

    # Calculate gap energy for the material
    delta = delta0(Tc)

    R2 = np.power(2*delta,2)/(2*N0*tau_0*np.power(c.k_B*Tc,3))
    
    nstar = np.power(R2*tau_max,-1).to(np.power(u.micron,-3))
    
    return nstar

#%%
filt_k,filt_t = np.loadtxt(module_path+'BandpassFilter.txt',skiprows=1,delimiter=',',unpack=True) # annoying, but gotta keep BandpassFilter.txt accessible
indx = np.argsort(filt_k) # sort the filter transmission spectrum by wavenumber
freq_filt = (filt_k[indx]*u.k).to(u.THz,equivalencies=u.spectral()) # apply the sort to the wavenumber points, convert from k=icm to THz
trans_filt = filt_t[indx] # apply the sort to the transmission points
from scipy.interpolate import CubicSpline 
redfilt = CubicSpline(freq_filt.value,trans_filt) # redfilt is now a function that takes in a freq in THz, returns a transmission value at that freq

''' Function to calculate the incident power emitted by a blackbody source at temperature T with a band-defining filter filt 
    and blocking filter with uniform transmission trans
     Inputs -- T_BB: blackbody temperature in K (can be an array); T_BB=0 returns P=0
               trans: transmission of the uniform blocking filter (unitless); default is 1 which corresponds to no filter, trans=0 returns P=0
               filt: spline function defining the transmission spectrum of the band-defining filter filt(nu_opt)=T, where nu_opt is in THz and T is unitless; default is the 350 micron QMC Spifi filter in Red cryostat
    Outputs -- P: power incident on the detectors in pW
'''
from astropy.modeling.blackbody import blackbody_nu
from scipy.integrate import trapz

def TBB_to_Pinc(T_BB,trans=1,filt=redfilt):
    # Treat scalars and arrays separately
    if np.size(T_BB)==1:
        if T_BB==0 or trans==0: P=0*u.pW # Make sure the scalar zero case returns zero power
        else: # Calculate BB power for the scalar case
            f = np.linspace(u.THz*filt.x.min(),u.THz*filt.x.max(),1000) # integrate over the range where we have filter transmission data
            B_nu = blackbody_nu(f,T_BB).to(u.W/u.m**2/u.steradian/u.Hz) # calculate the blackbody flux (per sr) at this blackbody temp for each frequency
            P_nu = B_nu*np.power(f.to(u.m,equivalencies=u.spectral()),2)*u.steradian  # Swapping units
            tr = trans*filt(f.to(u.THz))*u.dimensionless_unscaled # power passes through both filters before getting to the detectors
            P = trapz(P_nu*tr,x=f).to(u.pW) # Integrate the spectrum to get total power

    else:    
        P = np.zeros(len(T_BB))*u.pW # Return one power for each value of T_BB
        for ind,temp in enumerate(T_BB):
            if temp==0 or trans==0: P[ind] = 0*u.pW # again, treat the zero case
            else:
                f = np.linspace(freq_filt.min(),freq_filt.max(),1000) # integrate over the range where we have filter transmission data
                B_nu = blackbody_nu(f,temp).to(u.W/u.m**2/u.steradian/u.Hz) # calculate the blackbody flux (per sr) at this blackbody temp for each frequency
                P_nu = B_nu*np.power(f.to(u.m,equivalencies=u.spectral()),2)*u.steradian # Swapping units
                tr = trans*filt(f.to(u.THz))*u.dimensionless_unscaled # power passes through both filters before getting to the detectors
                P[ind] = trapz(P_nu*tr,x=f).to(u.pW) # Integrate the spectrum to get total power
    return P
#%%
''' Function to calculate the number density of thermal quasiparticles in the material at a stage temperature Tstage
     Inputs -- Tstage: Stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              N0: Density of states for the superconductor material (default is thin-film Al value)
    Outputs -- n_th: number density of thermal quasiparticles (in microns^-3)
'''
def nth(Tstage,Tc,N0=N0_Al):
    # Get our units in order
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc = Tc*u.K
    Tc = Tc.to(u.K)
    
    # Calculate gap energy for the material
    delta = delta0(Tc)
    
    # Calculate the thermal quasiparticle density
    n_th = (2*N0*(np.sqrt(2*np.pi*c.k_B*Tstage*delta)).to(u.eV)*np.exp(-delta/(c.k_B*Tstage))).to(np.power(u.micron,-3))
    
    return n_th

#%%
''' Function to calculate the number of total quasiparticles from thermal and photon generation/recombination
    Inputs -- Tstage: stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              T_BB: Blackbody temperature in K (can be an array)
              V: optically-active volume of the KID pixel in microns^3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_pb: pair breaking efficiency in the superconductor (unitless)
              eta_opt: optical efficiency *of the detector* (unitless)
              trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
    Output -- n_qp: number density of quasiparticles in the superconductor
'''
def nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans=1):
    # Get out units in order
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc*=u.K
    Tc = Tc.to(u.K)
    
    if not hasattr(T_BB,'unit'): T_BB*=u.K
    T_BB = T_BB.to(u.K)
    
    if not hasattr(V,'unit'): V = V*np.power(u.micron,3)
    V = V.to(np.power(u.micron,3))
    
    if not hasattr(tau_max,'unit') or tau_max.unit==u.dimensionless_unscaled: tau_max*=u.microsecond
    tau_max = tau_max.to(u.microsecond)
    
    if not hasattr(n_star,'unit') or n_star.unit==u.dimensionless_unscaled: n_star*=np.power(u.micron,-3)
    n_star = n_star.to(np.power(u.micron,-3))
    
    eta_pb*=u.dimensionless_unscaled

    trans*=u.dimensionless_unscaled
    
    # Calculate gap energy for the material
    delta = delta0(Tc)
    
    # Calculate the number of thermal quasiparticles
    n_th = nth(Tstage,Tc)
    
    # Calculate the power incident on the detector
    P_inc = TBB_to_Pinc(T_BB,trans)
    
    # Calculate the power absorbed by the detector
    P_abs = eta_opt*P_inc
    
    # Calculate the total quasiparticle density in the detector
    n_qp = (-n_star + ((n_star+n_th)**2 + (2*n_star*eta_pb*P_abs*tau_max)/(delta*V))**0.5).to(np.power(u.micron,-3))
    
    return n_qp

#%%
''' Function to calculate the electron temperature in a superconductor at a given stage temperature Tstage and under a given blackbody load
    The math to get to the analytical form of the equation comes from setting nth(T)=nqp(Pabs,T) and solving for T
    Inputs -- Tstage: Stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              T_BB: Blackbody temperature in K (can be an array)
              V: optically-active volume of the KID pixel in microns^3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_pb: pair breaking efficiency in the superconductor (unitless)
              eta_opt: optical efficiency *of the detector* (unitless)
              trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
              N0: Density of states for the superconductor material (default is thin-film Al value)
    Output -- T_elec: electron temperature in K
'''
from scipy.special import lambertw
def Telec(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans=1,N0=N0_Al):
    # Calculate the gap energy in the material
    delta = delta0(Tc)
    
    # Calculate the total number density of quasiparticles in the material
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # M is just a convenient parameterization 
    M = np.power(2*N0*delta,-1) * np.power(2*np.pi,-0.5) * n_qp
    M = M.to(u.dimensionless_unscaled)
    arg = (2/np.power(M,2)).value
    
    # The analytical solution to the inverse function of nth(T)=nqp(Pabs,T)
    y = 0.5*lambertw(arg)
    
    # Unpack the parameterized variable
    T_elec = delta/(c.k_B*y)
    T_elec = np.abs(T_elec).to(u.K)
    
    return T_elec

#%%
''' Function to calculate S1. 
    Inputs -- Tstage: Stage temperature in K (can be an array)
              f: KID pixel *Resonance Frequency* in MHz or similar
              Tc: Critical temperature of the superconductor in K
    Outputs -- S_1: (unitless)
'''
from scipy.special import kv # Bessel function
def S1(Tstage,f,Tc):
    # Get our units in order
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

    # (Just for convenience)
    beta = 1/(c.k_B*Tstage)
    
    # Calculate the gap energy in the material
    delta = delta0(Tc)
    
    # (Again, for convenience; note that sinh and kv don't play well with astropy quantities)
    arg = (0.5*beta*c.h*f).to(u.dimensionless_unscaled).value
    
    # Calculate the S_1 parameter
    S_1 = (2/np.pi)*np.sqrt((2*delta*beta)/(np.pi))*np.sinh(arg)*kv(0,arg)
    S_1 = S_1.to(u.dimensionless_unscaled)
    
    return S_1

#%%
''' Function to calculate S2. 
    Inputs -- Tstage: Stage temperature in K (can be an array)
              f: KID pixel *Resonance Frequency* in MHz or similar
              Tc: Critical temperature of the superconductor in K
    Outputs -- S_2: (unitless)
'''
from scipy.special import iv # Bessel function
def S2(Tstage,f,Tc):
    # Get our units in order 
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

    # (Just for convenience)
    beta = 1/(c.k_B*Tstage)
    
    # Calculate the gap energy in the material
    delta = delta0(Tc)
    
    # (Again, for convenience; note that iv doesn't play well with astropy quantities)
    arg = (0.5*beta*c.h*f).to(u.dimensionless_unscaled).value

    # Calculate the S_2 parameter
    S_2 = 1 + np.sqrt((2*delta*beta)/(np.pi))*np.exp(-arg)*iv(0,arg)
    S_2 = S_2.to(u.dimensionless_unscaled)
    
    return S_2

#%%
''' Function to calculate the quasiparticle lifetime in the superconductor given the quasiparticle number density 
    and the constants n_star and tau_max for the material
    Inputs -- n_qp: quasiparticle number density in microns^-3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds 
    Output -- tau_qp: quasiparticle lifetime in the material in microseconds
'''
def tauqp(n_qp,n_star,tau_max):
    # Get our units in order
    if not hasattr(n_qp,'unit'): n_qp*=np.power(u.micron,-3)
    n_qp = n_qp.to(np.power(u.micron,-3))
    
    if not hasattr(tau_max,'unit'): tau_max*=u.microsecond
    tau_max = tau_max.to(u.microsecond)
    
    if not hasattr(n_star,'unit'): n_star*=np.power(u.micron,-3)
    n_star = n_star.to(np.power(u.micron,-3))
    
    # Calculate the quasiparticle lifetime
    tau_qp = (tau_max/(1+n_qp/n_star)).to(u.microsecond)
    
    return tau_qp

#%%
''' Function to calculate the generation rate of thermal quasiparticles
    Inputs -- Tstage: Stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              V: optically-active volume of the KID pixel in microns^3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds 
    Output -- gamma_th: thermal quasiparticle generation rate in microseconds^-1
'''
def gammath(Tstage,Tc,V,n_star,tau_max):
    # Get our units in order
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc*=u.K
    Tc = Tc.to(u.K)
    
    if not hasattr(V,'unit'): V = V*np.power(u.micron,3)
    V = V.to(np.power(u.micron,3))

    if not hasattr(n_star,'unit'): n_star*=np.power(u.micron,-3)
    n_star = n_star.to(np.power(u.micron,-3))
    
    if not hasattr(tau_max,'unit'): tau_max*=u.microsecond
    tau_max = tau_max.to(u.microsecond)
    
    # Calculate the number density of thermal quasiparticles
    n_th = nth(Tstage,Tc)
    
    # Calculate the thermal quasiparticle generation rate
    gamma_th = (((n_th*V)/2)*(np.power(tau_max,-1)+np.power(tauqp(n_th,n_star,tau_max),-1))).to(np.power(u.microsecond,-1))
    
    return gamma_th

#%%
''' Function to calculate the recombination rate of quasiparticles in a superconductor
    Inputs -- Tstage: Stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              T_BB: Blackbody temperature in K (can be an array)
              V: optically-active volume of the KID pixel in microns^3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_pb: pair breaking efficiency in the superconductor (unitless)
              eta_opt: optical efficiency *of the detector* (unitless)
              trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
    Output -- gamma_r: quasiparticle recombination rate in microseconds^-1
'''
def gammar(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans=1):
    # Get our units in order
    if not hasattr(Tstage,'unit'):
        if Tstage.any() > 1.5: Tstage = Tstage*u.mK
        else: Tstage = Tstage*u.K
    Tstage = Tstage.to(u.K)
                       
    if not hasattr(Tc,'unit'): Tc*=u.K
    Tc = Tc.to(u.K)
    
    if not hasattr(T_BB,'unit'): T_BB*=u.K
    T_BB = T_BB.to(u.K)
    
    if not hasattr(V,'unit'): V = V*np.power(u.micron,3)
    V = V.to(np.power(u.micron,3))

    if not hasattr(tau_max,'unit'): tau_max*=u.microsecond
    tau_max = tau_max.to(u.microsecond)
    
    if not hasattr(n_star,'unit'): n_star*=np.power(u.micron,-3)
    n_star = n_star.to(np.power(u.micron,-3))

    # Calculate the total number density of quasiparticles
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Calculate the quasiparticle lifetime
    tau_qp = tauqp(n_qp,n_star,tau_max)
    
    # Calculate the quasiparticle recombination rate
    gamma_r = (((n_qp*V)/2)*(np.power(tau_max,-1)+np.power(tau_qp,-1))).to(np.power(u.microsecond,-1))
    return gamma_r

#%%
''' Function to calculate the fractional frequency shift relative to the zero temperature & loading state
    of a resonator under a given set of thermal & optical conditions
    Inputs -- alpha: kinetic inductance fraction of the inductor (unitless)
              f: KID pixel *Resonance Frequency* in MHz or similar
              Tstage: Stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              T_BB: Blackbody temperature in K (can be an array)
              V: optically-active volume of the KID pixel in microns^3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_pb: pair breaking efficiency in the superconductor (unitless)
              eta_opt: optical efficiency *of the detector* (unitless)
              trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
              N0: Density of states for the superconductor material (default is thin-film Al value)
    Output -- x_MB: fractional freq shift relative to the T = nqp = 0 state (unitless)
'''

def xMB(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans=1,N0=N0_Al):
    # Calculate the gap energy
    delta = delta0(Tc)
    
    # Calculate the number density of quasiparticles
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Calculate the electron temperature in the material
    T_elec = Telec(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans,N0)
    
    # Calculate the fractional frequency shift, using T_elec to calculate S2 instead of the physical temperature
    x_MB = -n_qp*(alpha*S2(T_elec,f,Tc))/(4*N0*delta)
    x_MB = x_MB.to(u.dimensionless_unscaled)
    
    return x_MB

#%%
''' Function to calculate the responsivity*optical efficiency
    of a resonator under a given set of thermal & optical conditions
    Inputs -- alpha: kinetic inductance fraction of the inductor (unitless)
              f: KID pixel *Resonance Frequency* in MHz or similar
              Tstage: Stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              T_BB: Blackbody temperature in K (can be an array)
              V: optically-active volume of the KID pixel in microns^3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_pb: pair breaking efficiency in the superconductor (unitless)
              eta_opt: optical efficiency *of the detector* (unitless)
              trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
              N0: Density of states for the superconductor material (default is thin-film Al value)
    Output -- x_MB: fractional freq shift relative to the T = nqp = 0 state (unitless)
'''

def dxdPabs(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans=1,N0=N0_Al):
    # Calculate the gap energy
    delta = delta0(Tc)
    
    # Calculate the number density of quasiparticles
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Calculate the quasiparticle lifetime in the material
    tau_qp = tauqp(n_qp,n_star,tau_max)
    
    # Calculate the electron temperature in the material
    T_elec = Telec(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans,N0)
    
    # Calculate the fractional frequency shift, using T_elec to calculate S2 instead of the physical temperature
    dx_dPabs = eta_opt*((alpha*S2(T_elec,f,Tc))/(4*N0*delta))*((eta_pb*tau_qp)/(delta*V))
    dx_dPabs = dx_dPabs.to(np.power(u.W,-1))
    
    return dx_dPabs
    
#%%
''' Function to calculate 1/Q_MB of a resonator under a given set of thermal & optical conditions
    Inputs -- alpha: kinetic inductance fraction of the inductor (unitless)
              f: KID pixel *Resonance Frequency* in MHz or similar
              Tstage: Stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              T_BB: Blackbody temperature in K (can be an array)
              V: optically-active volume of the KID pixel in microns^3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_pb: pair breaking efficiency in the superconductor (unitless)
              eta_opt: optical efficiency *of the detector* (unitless)
              trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
              N0: Density of states for the superconductor material (default is thin-film Al value)
    Output -- Q_inv_MB: 1/Q_MB (unitless)

'''
def QinvMB(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans=1,N0=N0_Al):
    # Calculate the gap energy
    delta = delta0(Tc)
    
    # Calculate the total quasiparticle density
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Calculate the electron temperature
    T_elec = Telec(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans,N0)

    # Calculate S1 using the electron temperature rather than the physical temperature
    S_1 = S1(T_elec,f,Tc)
    
    # Calculate the inverse resonator quality factor
    Q_inv_MB = n_qp*(alpha*S_1)/(2*N0*delta)
    Q_inv_MB = Q_inv_MB.to(u.dimensionless_unscaled)
    
    return Q_inv_MB


#%%
''' Function to calculate fractional frequency (white) noise Sxx under a given set of temperature and loading conditions
    Inputs -- alpha: kinetic inductance fraction of the inductor (unitless)
              f: KID pixel *Resonance Frequency* in MHz or similar
              Tstage: Stage temperature in K (can be an array)
              Tc: Critical temperature of the superconductor in K
              T_BB: Blackbody temperature in K (can be an array)
              V: optically-active volume of the KID pixel in microns^3
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_pb: pair breaking efficiency in the superconductor (unitless)
              nu_opt: frequency of optical photons (in Hz or similar)
              eta_opt: optical efficiency *of the detector* (unitless)
              trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
              N0: Density of states for the superconductor material (default is thin-film Al value)
    Output -- S_xx: fractional frequency noise (in Hz^-1)

'''

def Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt,trans=1,N0=N0_Al):
    # Calculate the gap energy 
    delta = delta0(Tc)

    # Calculate the incident and absorbed powers
    P_inc = TBB_to_Pinc(T_BB,trans)
    P_abs = eta_opt*P_inc

    # Calculate the effective electron temperature
    T_elec = Telec(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans,N0)
    
    # Calculate S2 using the effective electron temperature rather than the physical temperature
    S_2 = S2(T_elec,f,Tc)
    
    # Calculate the photon occupation number
    n_gamma = ngamma(nu_opt,T_BB)
    
    # Calculate the total quasiparticle number density
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Calculate the quasiparticle lifetime in the material
    tau_qp = tauqp(n_qp,n_star,tau_max)
    
    # Calculate the thermal quasiparticle generation rate 
    gamma_th = gammath(Tstage,Tc,V,n_star,tau_max)
    
    # Calculate the quasiparticle recombination rate
    gamma_r = gammar(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Phew. Now calculate the actual S_xx white noise level
    S_xx = np.power(alpha*S_2/(4*N0*delta),2) * (np.power(eta_pb*tau_qp/(delta*V),2)*(2*c.h*nu_opt*P_abs)*(1+n_gamma) + 4*np.power(tau_qp,2)*(gamma_th+gamma_r)/np.power(V,2))
    S_xx = S_xx.to(np.power(u.Hz,-1))
    
    return S_xx

#%%
def Sxx_G_photon(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt,trans=1,N0=N0_Al):
    # Calculate the gap energy 
    delta = delta0(Tc)

    # Calculate the incident and absorbed powers
    P_inc = TBB_to_Pinc(T_BB,trans)
    P_abs = eta_opt*P_inc

    # Calculate the effective electron temperature
    T_elec = Telec(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans,N0)
    
    # Calculate S2 using the effective electron temperature rather than the physical temperature
    S_2 = S2(T_elec,f,Tc)
    
    # Calculate the photon occupation number
    n_gamma = ngamma(nu_opt,T_BB)
    
    # Calculate the total quasiparticle number density
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Calculate the quasiparticle lifetime in the material
    tau_qp = tauqp(n_qp,n_star,tau_max)

    S_xx_G_photon = np.power(alpha*S_2/(4*N0*delta),2) * np.power(eta_pb*tau_qp/(delta*V),2) * (2*c.h*nu_opt*P_abs)*(1+n_gamma)

    return S_xx_G_photon.to(np.power(u.Hz,-1))

#%%
def Sxx_R_photon(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt,trans=1,N0=N0_Al):
    # Calculate the gap energy 
    delta = delta0(Tc)

    # Calculate the incident and absorbed powers
    P_inc = TBB_to_Pinc(T_BB,trans)
    P_abs = eta_opt*P_inc

    # Calculate the effective electron temperature
    T_elec = Telec(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans,N0)
    
    # Calculate S2 using the effective electron temperature rather than the physical temperature
    S_2 = S2(T_elec,f,Tc)
    
    # Calculate the photon occupation number
    n_gamma = ngamma(nu_opt,T_BB)
    
    # Calculate the total quasiparticle number density
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Calculate the quasiparticle lifetime in the material
    tau_qp = tauqp(n_qp,n_star,tau_max)

    S_xx_R_photon = np.power(alpha*S_2/(4*N0*delta),2) * np.power(eta_pb*tau_qp/(delta*V),2) * 4*P_abs*delta/eta_pb

    return S_xx_R_photon.to(np.power(u.Hz,-1))
    
#%%
def Sxx_GR_th(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt,trans=1,N0=N0_Al):
    delta = delta0(Tc)

    # Calculate the effective electron temperature
    T_elec = Telec(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans,N0)
    
    # Calculate S2 using the effective electron temperature rather than the physical temperature
    S_2 = S2(T_elec,f,Tc)
        
    # Calculate the total quasiparticle number density
    n_qp = nqp(Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Calculate the quasiparticle lifetime in the material
    tau_qp = tauqp(n_qp,n_star,tau_max)

    # Calculate the thermal quasiparticle generation rate 
    gamma_th = gammath(Tstage,Tc,V,n_star,tau_max)
    
    S_xx_GR_th = np.power(alpha*S_2/(4*N0*delta),2) * np.power(2*tau_qp/V,2) * gamma_th
    
    return S_xx_GR_th.to(np.power(u.Hz,-1))
    
    