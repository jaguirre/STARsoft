#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 10:07:23 2017

@author: tashaleebillings
"""

"""
 NAME:
     Line Intensity Mapping CII Regions

 PURPOSE:
     To predict the clustering and evolution behavior of galaxies using the mean intensity with the Far-IR line intensity of CII regions 
     The model for the signal has just a couple of components:
         - How the luminosity functions evolve with redshift (Equation 4, and the parameters that follow it, plus Equation A1)
         - How to map the line luminosity onto the L_IR (Equation 6)
         - The matter power spectrum (derived from CAMB (see the beginning of Sec 2.1)
         - The cosmological parameters y, D_A, and D_L (can be calculated from astropy.cosmology)
         - Then the expected power spectrum is just the sum of Equations 1 and 3.

    PAPER USED: http://adsabs.harvard.edu/abs/2014ApJ...793..116U
    
 CATEGORY:

 CALLING SEQUENCE:

 INPUTS:
     z -- redshift  
     L_IR -- Luminosity of an Infrared Galaxy (IR Luminosity) [W] (Local galaxies have luminosities between 1e8*L_0 and 1e13*L_0)
       
 OPTIONAL INPUTS:

 KEYWORD PARAMETERS:
     B11Lum -- Luminosity Function as a function of L_IR and z [W]
     Lstar -- Luminosity as function of redshit [W]
     Phistar -- Luminosity Function of any type of galaxy [1/dex/Mpc^3] 
         dex(an order of magnetude) --> 10e0.1 = 0.1 dex
     L_0 -- Luminosity of Sun [W] {default/actual value is 3.846e26}
     rPhi -- redshift dependent
     rL -- redshift dependent
     beta = faint-end power law slope of the luminosity function {default value is 1.223}
     xi = bright-end Gaussian width of the luminosity function {default values is 0.406}

 OUTPUTS:
     ['Lstar']
     ['Phistar']
     ['B11Lum']
     
 RESTRICTIONS:

 EXAMPLE:

 MODIFICATION HISTORY:  
       


"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sf
from scipy.integrate import quad
import astropy.units as u
import astropy.constants as c
from astropy.cosmology import Planck15
import cosmolopy.perturbation as cp
"""
Constants
"""
zbreak1 = 0.879 # This break in the redshift is a fitted parameter
zbreak2 = 2 # This break in the redshift is a fixed value
L_0 = c.L_sun
beta = 1.223
xi = 0.406

def phistar(z):
    if z < zbreak1:
        rPhi = 0.774
    elif zbreak1 < z and z < zbreak2:
        rPhi = -6.246
    else:
        rPhi = -0.919
#    Mpc = 1e6*u.pc   
    answer = 3.234e-3*(1/u.dex/(np.power(u.Mpc,3)))*np.power(1+z,rPhi)
    
    return answer

def lstar(z):
    if z < zbreak1:
        rL = 2.931
    elif zbreak1 < z and z < zbreak2:
        rL = 4.737
    else:
        rL = 0.145
    
    answer = 2.377e10*np.power(1+z,rL)
    
    return answer

#%%    
def b11_IRLum_Func(L_IR,z): #Bethermin 11 Function / Schechter Luminosity Function
    Phistar = phistar(z).value #Normalization Point
    Lstar = lstar(z)
    #ll = np.linspace(1e8,1e13, num=10000)
    #L_IR = ll
    B11Lum = Phistar*np.power(L_IR/Lstar,1-beta)*np.exp(-0.5/(np.power(xi,2))*np.log(1+(L_IR/Lstar)))
    
    answer = {'Lstar':Lstar, 'Phistar':Phistar, 'B11Lum':B11Lum}
    
    return answer

def Si_integrand(L_IR,z): #eq7 B.D. Uzgil
    lamb_rest = (158*u.micron.to(u.km))*u.km #[mincron]
    A,B,sigma_A,sigma_B = 0.89,2.44,0.03,0.07
    D_Aco = Planck15.angular_diameter_distance(z).value #Angular Diameter Distance: Gives info on the physical and angular size of a distant object [Mpc]
    D_L = np.power(1+z,2)*D_Aco # Luminosity Distance
# The ratio of D_Aco and D_L is the square inverse (1+z)
    
    H = Planck15.H0.value*np.sqrt(Planck15.Om(z)*np.power(1+z,3) +Planck15.Ok(z)*np.power(1+z,2) +Planck15.Ode(z)) # Hubble Parameter as a function of redshift
    y = lamb_rest*(np.power(1+z,2)/H) #[Mpc * s]
    x_plus = A+sigma_A*np.log10(L_IR)-(B+sigma_B) # Eq 5&6 B.D. Uzgil
    x_minus = A-sigma_A*np.log10(L_IR)-(B-sigma_B) # Eq 5&6 B.D. Uzgil
    f_iplus = (np.power(10,x_plus)/L_IR) 
    f_iminus = (np.power(10,x_minus)/L_IR)

    
    func_plus = b11_IRLum_Func(L_IR,z)['B11Lum']*((f_iplus*L_IR)/(4*np.pi*np.power(D_L,2)))*y*np.power(D_Aco,2)
    func_minus = b11_IRLum_Func(L_IR,z)['B11Lum']*((f_iminus*L_IR)/(4*np.pi*np.power(D_L,2)))*y*np.power(D_Aco,2)

    answer = {'plus':func_plus, 'minus':func_minus}

    return answer

def P_shot_integrand(L_IR,z):
    lamb_rest = 158 #[mincron]
    A,B,sigma_A,sigma_B = 0.89,2.44,0.03,0.07
    D_Aco = Planck15.angular_diameter_distance(z).value #Angular Diameter Distance: Gives info on the physical and angular size of a distant object [Mpc]
    D_L = np.power(1+z,2)*D_Aco # Luminosity Distance
# The ratio of D_Aco and D_L is the square inverse (1+z)
    
    H = Planck15.H0.value*np.sqrt(Planck15.Om(z)*np.power(1+z,3) +Planck15.Ok(z)*np.power(1+z,2) +Planck15.Ode(z)) # Hubble Parameter as a function of redshift
    y = lamb_rest*(np.power(1+z,2)/H)
    x_plus = A+sigma_A*np.log10(L_IR)-(B+sigma_B) # Eq 5&6 B.D. Uzgil
    x_minus = A-sigma_A*np.log10(L_IR)-(B-sigma_B) # Eq 5&6 B.D. Uzgil
    f_iplus = (np.power(10,x_plus)/L_IR)
    f_iminus = (np.power(10,x_minus)/L_IR)

    
    func_plus = b11_IRLum_Func(L_IR,z)['B11Lum']*np.power(((f_iplus*L_IR)/(4*np.pi*np.power(D_L,2)))*y*np.power(D_Aco,2),2)
    func_minus = b11_IRLum_Func(L_IR,z)['B11Lum']*np.power(((f_iminus*L_IR)/(4*np.pi*np.power(D_L,2)))*y*np.power(D_Aco,2),2)

    answer = {'plus':func_plus, 'minus':func_minus}

    return answer

def P_cluster(L_IR,k,z):
    b = [2,2.3,2.6,2.9]# The linear bias is normally a function of z but in the paper 4 values are given
    cosmology = {'omega_M_0':Planck15.Om0, 'omega_lambda_0':Planck15.Ode0, 'omega_k_0':Planck15.Ok0, 
                 'h':Planck15.H0.value/100, 'omega_b_0' : Planck15.Ob0,'omega_n_0' : 0.0,'n' : 1.0,
                 'sigma_8' : 0.9,'N_nu' : 0,'baryonic_effects' : False}
    P_deltadelta = cp.power_spectrum(k,0,**cosmology) # This is probably the linear 
    
    power_plus = np.power(Si_integrand(L_IR,z)['plus']*b,2)*P_deltadelta
    power_minus = np.power(Si_integrand(L_IR,z)['minus']*b,2)*P_deltadelta

    answer = {'plus':power_plus,'minus':power_minus}
    
    return answer

#Planck15.H
#%%
#NOT GETTING GREAT STARTING POINTS BUT SHAPE LOOKS GOOD. Z=1 IS AN ISSUE
#k = np.logspace(-3,0)
#z=0
#Phistar = phistar(z).value
#Lstar = lstar(z)
#ll = np.linspace(1e8,1e13, num=10000)
#L_IR = ll
#B11Lum = Phistar*np.power(L_IR/Lstar,1-beta)*np.exp(-0.5/(np.power(xi,2))*np.log(1+(L_IR/Lstar)))
#plt.loglog(ll,B11Lum,'-k',label='z= '+str(z))
#plt.legend()

#%%
#np.trapz(,x=)


#cosmology = {'omega_M_0':Planck15.Om0, 'omega_lambda_0':Planck15.Ode0, 'omega_k_0':Planck15.Ok0, 'h':Planck15.H0.value/1e2, 'omega_b_0' : Planck15.Ob0,'omega_n_0' : 0.0,'n' : 1.0,'sigma_8' : 0.9,'N_nu' : 0,'baryonic_effects' : False}
#cp.norm_power(**cosmology)
l_ir = np.logspace(8,13,num=10000)
plt.figure('Bethermin 11 Function')
plt.clf()
plt.loglog(l_ir,b11_IRLum_Func(l_ir,z=0)['B11Lum'],'k', label = 'z=0')
plt.loglog(l_ir,b11_IRLum_Func(l_ir,z=1)['B11Lum'],'r', label = 'z=1')
plt.loglog(l_ir,b11_IRLum_Func(l_ir,z=1.5)['B11Lum'],'b', label = 'z=1.5')
plt.loglog(l_ir,b11_IRLum_Func(l_ir,z=2)['B11Lum'],'g', label = 'z=2')
plt.loglog(l_ir,b11_IRLum_Func(l_ir,z=3)['B11Lum'],'m', label = 'z=3')

plt.ylabel(r'$\Phi$')
plt.xlabel('L_IR')
plt.legend()#(handles,labels,bbox_to_anchor=(.95,.6),borderaxespad=0.)
plt.show()








