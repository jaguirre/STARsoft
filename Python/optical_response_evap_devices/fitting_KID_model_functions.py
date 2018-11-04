# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 20:15:08 2018
@author: Alyssa
"""

import numpy as np
from astropy import units as u
from astropy import constants as c

import KID_model_functions as kids

#%%
''' Function specifically for use to fit x vs Tstage dark data for the free parameters alpha, Tc, and dx
    Inputs -- data: data[0] is an array of stage temperatures in K
                    data[1] is an individual KID resonance frequency in MHz or similar
             alpha: kinetic inductance fraction of the inductor (unitless)
             Tc: Critical temperature of the superconductor in K
             df: nuisance parameter to account for the fact that we can't measure at zero temperature and optical loading. 
                 *** df is in units inverse to f, but needs to be given as unitless for the fitter ***
                 x_tot = x_MB + df/f = x_MB + dx -> dx = f*df
    Output -- x_tot: total fractional freq shift (unitless)
'''
def x_dark_fit(data,alpha,Tc,df):
    # Unpack the data input
    Tstage = data[0]
    f = data[1]
    
    # Calculate x_MB for the given parameters
    x_MB = kids.xMB(alpha,f,Tstage,Tc,T_BB=1,V=1,n_star=1,tau_max=1,eta_pb=1,eta_opt=0,trans=0) 
    
    # Correct for the fact that scipy.io.curve_fit doesn't deal well with units
    df = df*np.power(f.unit,-1)
    dx = (df*f).to(u.dimensionless_unscaled)
    
    # add in the nuisance parameter dx = df*f
    x_tot = x_MB + dx
    return x_tot

#%%
''' Function specifically for use to fit 1/Q vs Tstage dark data for the free parameters alpha, Tc, and Q_inv0
    Inputs -- data: data[0] is an array of stage temperatures in K
                    data[1] is an individual KID resonance frequency in MHz or similar
             alpha: kinetic inductance fraction of the inductor (unitless)
             Tc: Critical temperature of the superconductor in K
             Q_inv0: nuisance parameter to account for the coupling and internal quality factors
                     Q_inv = Q_inv_MB + Q_inv0 -> 1/Q = 1/Q_MB + 1/Q_0
    Output -- Q_inv: Resonator quality factor (unitless)
'''

def Qinv_dark_fit(data,alpha,Tc,Q_inv0=0):
    # Unpack the input parameter
    Tstage = data[0]
    f = data[1]
    Q_inv0*=u.dimensionless_unscaled
    
    # Calculate the Q_MB value for the given conditions
    Q_inv_MB = kids.QinvMB(alpha,f,Tstage,Tc,T_BB=1,V=1,n_star=1,tau_max=1,eta_pb=1,eta_opt=0,trans=0)
    
    # Add in the Q_inv0 shift
    Q_inv = Q_inv_MB + Q_inv0
    return Q_inv

#%%
''' Dummy function to simultaneously fit dark x vs Tstage and Qinv vs Tstage data
    Inputs -- data: data[0] is an array of stage temperatures in K
                    data[1] is an individual KID resonance frequency in MHz or similar
             alpha: kinetic inductance fraction of the inductor (unitless)
             Tc: Critical temperature of the superconductor in K
             df: nuisance parameter to account for the fact that we can't measure at zero temperature and optical loading. 
             Q_inv0: nuisance parameter to account for the coupling and internal quality factors
    Output -- glom: concatenation of x and Qinv values (should always be of equal length)
'''
def x_Qinv_dark_simulfit(data,alpha,Tc,df,Q_inv0):
    # Calculate x and Qinv using their fit functions
    x = x_dark_fit(data,alpha,Tc,df)
    Qinv = Qinv_dark_fit(data,alpha,Tc,Q_inv0)

    # Mash them together so the fitter can return a 1D vector
    glom = np.concatenate((x,Qinv))
    return glom

#%%
''' Function specifically to fit x vs T_BB(->P_inc)
    Inputs -- data: data[0] = T_BB: Blackbody temperature in K (should be an array)
                    data[1] = alpha: kinetic inductance fraction of the inductor (unitless)
                    data[2] = f: KID pixel *Resonance Frequency* in MHz or similar
                    data[3] = Tstage: Stage temperature in K (should be scalar)
                    data[4] = Tc: Critical temperature of the superconductor in K
                    data[5] = V: optically-active volume of the KID pixel in microns^3
                    data[6] = eta_pb: pair breaking efficiency in the superconductor (unitless)
                    data[7] = nu_opt: frequency of optical photons (in Hz or similar)
                    data[8] = trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
                    data[9] =  N0: Density of states for the superconductor material (default is thin-film Al value)              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_opt: optical efficiency *of the detector* (unitless)
              df: nuisance parameter to account for the fact that we can't measure at zero temperature and optical loading. 
                 *** df is in units inverse to f, but needs to be given as unitless for the fitter ***
                 x_tot = x_MB + df/f = x_MB + dx -> dx = f*df
    Output -- x_opt: fractional frequency shift relative to zero temp and loading
'''
def x_opt_fit(data,n_star,tau_max,eta_opt,df):
    # Unpack the input parameter
    T_BB,alpha,f,Tstage,Tc,V,eta_pb,nu_opt,trans,N0 = data
    
    # Correct for the fact that scipy.io.curve_fit doesn't deal well with units
    df = df*np.power(f.unit,-1)
    dx = (df*f).to(u.dimensionless_unscaled)

    # Calculate x_opt
    x_MB = kids.xMB(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt,trans)
    
    # Add in the nuisance parameter dx
    x_opt = (x_MB + dx).to(u.dimensionless_unscaled)
    
    return x_opt

#%%
''' Function specifically to fit Sxx vs T_BB(->P_inc)
    Inputs -- data: data[0] = T_BB: Blackbody temperature in K (should be an array)
                    data[1] = alpha: kinetic inductance fraction of the inductor (unitless)
                    data[2] = f: KID pixel *Resonance Frequency* in MHz or similar
                    data[3] = Tstage: Stage temperature in K (should be scalar)
                    data[4] = Tc: Critical temperature of the superconductor in K
                    data[5] = V: optically-active volume of the KID pixel in microns^3
                    data[6] = eta_pb: pair breaking efficiency in the superconductor (unitless)
                    data[7] = nu_opt: frequency of optical photons (in Hz or similar)
                    data[8] = trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
                    data[9] =  N0: Density of states for the superconductor material (default is thin-film Al value)
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_opt: optical efficiency *of the detector* (unitless)
              Sxx_0: Extra constant to account for noise from unknown sources in Hz^-1 but unitless because scipy.optimize.curve_fit doesn't like units
    Output -- S_xx: fractional frequency noise (in Hz^-1 but unitless because scipy.optimize.curve_fit doesn't like units)
'''

def Sxx_fit(data,n_star,tau_max,eta_opt,Sxx_0):
    # Unpack the input parameter
    T_BB,alpha,f,Tstage,Tc,V,eta_pb,nu_opt,trans,N0 = data
    
    # Calculate Sxx
    S_xx = kids.Sxx(alpha,f,Tstage,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt,trans,N0)
    S_xx_tot = (S_xx + Sxx_0*(np.power(u.Hz,-1))).to(np.power(u.Hz,-1))
    S_xx_tot = S_xx_tot.value
    
    return S_xx_tot

#%%
''' Dummy function to simultaneously fit x_opt and Sxx as a function of TBB
    Inputs -- data: data[0] = T_BB: Blackbody temperature in K (should be an array)
                    data[1] = alpha: kinetic inductance fraction of the inductor (unitless)
                    data[2] = f: KID pixel *Resonance Frequency* in MHz or similar
                    data[3] = Tstage: Stage temperature in K (should be scalar)
                    data[4] = Tc: Critical temperature of the superconductor in K
                    data[5] = V: optically-active volume of the KID pixel in microns^3
                    data[6] = eta_pb: pair breaking efficiency in the superconductor (unitless)
                    data[7] = nu_opt: frequency of optical photons (in Hz or similar)
                    data[8] = trans: transmission of a blocking filter (unitless); trans=0 gives nqp=nth
                    data[9] =  N0: Density of states for the superconductor material (default is thin-film Al value)
              n_star: quasiparticle number constant in microns^-3
              tau_max: quasiparticle lifetime constant in microseconds
              eta_opt: optical efficiency *of the detector* (unitless)
              df: nuisance parameter to account for the fact that we can't measure at zero temperature and optical loading. 
                 *** df is in units inverse to f, but needs to be given as unitless for the fitter ***
                 x_tot = x_MB + df/f = x_MB + dx -> dx = f*df
              Sxx_0: Extra constant to account for noise from unknown sources in Hz^-1 but unitless because scipy.optimize.curve_fit doesn't like unit
    Output -- glom: concatenation of x_opt and S_xx values (should always be of equal length)
'''
def x_Sxx_opt_simulfit(data,n_star,tau_max,eta_opt,df,Sxx_0):
    # Calculate x_opt and Sxx separately using their fit functions
    x_opt = x_opt_fit(data[0],n_star,tau_max,eta_opt,df)
    S_xx_tot = Sxx_fit(data[1],n_star,tau_max,eta_opt,Sxx_0)
    
    # Mash them together so the fitter can return a 1D vector
    glom = np.concatenate((x_opt,S_xx_tot))
    
    return glom