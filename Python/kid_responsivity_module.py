#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:56:01 2017

@author: tashaleebillings
"""

"""
 NAME:
       KID_RESPONSIVITY

 PURPOSE:
       Compute the responsivity of a KID at finite temperature, and under a specified optical load. A number of associated quantities are computed as well. These are:
        - nth -- number density of quasiparticles for a superconductor in thermal equilibrium at the specified temperature [micron^-3]
        - nqp -- number density of quasiparticles for a superconductor at the specified temperature and optical load [micron^-3]
        - tau_qp -- quasiparticle lifetime for a superconductor at the specified temperature and optical load [microsec]
        - s1 -- S_1(w)
        - s2 -- S_2(w)
        - xr -- df/f, referenced to the nqp = 0 state
        - Qi_inv -- Qi^-1
        - r_x -- df/f responsivity [W^-1]
        - r_qinv -- Qi^-1 responsivity [W^-1]
        - sxx_gr -- GR noise, excluding the photon generation noise [Hz^-1]
        - sxx_gr0 -- alternative, simpler form of sxx_gr often found in the literature [Hz^-1]
        - sxx_gamma -- photon generation noise [Hz^-1]
        
        See 'KID_derivations.pdf' for details.
       
 CATEGORY:

 CALLING SEQUENCE:

 INPUTS:
       temp -- T [K]
       pabs -- absorbed power [pW]
       
 OPTIONAL INPUTS:

 KEYWORD PARAMETERS:
       tc -- Tc [K] {default is 1}
       N0 -- single spin density of states at the Fermi level [micron^-3 eV^-1] {default is 1.72e10 for Al}
       nstar - n* [micron^-3] {default is 100}
       tau_max -- tau_max [microsec] {default is 100}
       eta_pb -- pair-breaking efficiency {default is 0.7}
       vol -- inductor volume [micron^3] {default is 1}
       fr -- resonator frequency [MHz] {default is 100}
       alpha_k -- kinetic inductance fraction {default is 1}
       gamma_t -- thin film parameter {default is 1}
       nu_opt -- optical frequency for photon noise computation [GHz] {default is 250}
       n0_gamma -- photon occupation number in the detector, for photon noise computation {default is 0}

 OUTPUTS:
       struct -- a structure with the following fields
       .nth
       .nqp
       .tau_qp
       .s1
       .s2
       .xr
       .Qi_inv
       .r_x
       .r_qinv
       .sxx_gr
       .sxx_gr0
       .sxx_gamma

 RESTRICTIONS:

 EXAMPLE:

 MODIFICATION HISTORY:  
       07/24/17 SHD

-
***********************************************************************
"""

import numpy as np
import scipy.special as sf
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c
#from KID_Lifetime_Response_noise import nqp_tot

#Define various constants
k = c.k_B.value #[J/K]
h_p = c.h.value #[Js]
ev_joule = 1.6022e-19 #[eV/J]

def kid_responsivity(temp, pabs, tc = 1.2,N0 = 1.72e10,nstar = 1e2,tau_max = 1e-3,eta_pb = 0.7,vol = 1,fr = 100,alpha_k = 1,gamma_t = 1,nu_opt = 250,n0_gamma = 0): # Constants are found in McCarrick right below Eq. 26
#   Compute the ratio of Delta_0/k_B [K]
    d0_kB = 1.764*tc
    delta = d0_kB* k
    nqp_min=400

#   Compute nth
    N0 = N0*np.power(u.micron,-3)/u.eV # per microns^3 per eV
    N0 = (N0.to(np.power(u.micron,-3)/u.J)).value # per microns^3 per Joule but it's just the value
    nth = 2.*N0*k*np.sqrt(2.*np.pi*temp*d0_kB)*np.exp(-1.*d0_kB/temp)
    #nth = 2 * N0 * np.sqrt(2. * np.pi * delta * k * temp) * np.exp(-1*delta / (k * temp)) + nqp_min #Eq. 14 McCarrick and Eq. 51 Mauskopf Review

#   Compute nqp
#   This expression has a term of the form [sqrt(1 + eps) - 1], where eps is small when nth and pabs are small.
#   When eps is small this term is not computed accurately, and here we explicitly linearize it.
    eps = 2.*nth/nstar + np.power((nth/nstar),2) + 2.*eta_pb*pabs*1.e-12*tau_max*1.e-6/(nstar*vol*d0_kB*k)
    term = np.sqrt(1. + eps) - 1.
    indx = np.where(eps < 1.e-8) # if there are no matches indx is forced to be count
    #count = len(indx[0])
    #if count != 0:
    term[indx] = 0.5*eps[indx]
    nqp = nstar*term 

#   Compute tau_qp
    tau_qp = tau_max/(1. + nqp/nstar) #Eq. 22 McCarrick

#   Compute S1 and S2
    xi = (h_p*fr*1.e6)/(2.*k*temp) # Mauskopf Review
    s1 = (2./np.pi)*np.sqrt(2.*d0_kB/(np.pi*temp))*np.sinh(xi)*sf.kv(0,xi) #Eq. 64 Mauskopf and Eq. 16 Jonas Review
    #s2 = 1. + np.sqrt(2.*d0_kB/(np.pi*temp))*np.exp(-1.*xi)*sf.iv(0,xi) #Eq. 65 Mauskopf and Eq. 5 Jonas Review
    s2 = ((np.pi*delta)/(h_p*fr*1e6))*(1-nqp/(2*N0*delta))*(1. + np.sqrt(2.*d0_kB/(np.pi*temp))*np.exp(-1.*xi)*sf.iv(0,xi))

    
#    Check to how these fucntions behave
#    res_fr= np.arange(75,165)
#    xi = ((h_p*res_fr*1.e6)/(2.*k*0.2)).value # Mauskopf Review
#    s1 = (2./np.pi)*np.sqrt(2.*d0_kB/(np.pi*0.2))*np.sinh(xi)*sf.kv(0,xi)
#    s2 = 1. + np.sqrt(2.*d0_kB/(np.pi*0.2))*np.exp(-1.*xi)*sf.iv(0,xi) 
#    plt.plot(res_fr,beta)

#   Compute xr and Qi_inv
#   Note that xr refers to the frequency shift from the nqp = 0 state
    xr = -1.*alpha_k*gamma_t*s2*nqp/(4.*N0*d0_kB*k) #Eq. 18 McCarrick
    Qi_inv = -1.*xr*2.*s1/s2 #Eq. 63 Mauskopf and Eq. 10/12 McCarrick

#   Compute the frequency and Qinv responsivity
    r_x = -1.*alpha_k*gamma_t*s2/(4.*N0*d0_kB*k)*eta_pb*tau_qp*1.e-6/(d0_kB*k*vol) #Below Eq. 5 Steve Hailey-Dunsheath and Eq. 26 McCarrick: The detector responsivity for absorbed power
    r_qinv = -1.*r_x*2.*s1/s2
    
#   Compute Sxx_gr and Sxx_gr0
    tau_th = tau_max/(1. + nth/nstar) #quasiparticle lifetime for a superconductor in thermal equilibrium at the specified temperature [microsec]
    gamma_th = nth*vol/2.*(1./tau_max + 1./tau_th)*1.e6 #quasiparticle generation rate due to thermal fluctuations ;[sec^-1]
    gamma_r = nqp*vol/2.*(1./tau_max + 1./tau_qp)*1.e6 #quasiparticle recombination rate ;[sec^-1]
    sxx_gr = np.power((alpha_k*gamma_t*s2/(4.*N0*d0_kB*k)),2)*4.*np.power((tau_qp*1.e-6)/vol,2)*(gamma_th + gamma_r)
    sxx_gr0 = np.power((alpha_k*gamma_t*s2/(4.*N0*d0_kB*k)),2)*4.*nqp/vol*tau_qp*1.e-6

#   Compute Sxx_gamma
    sxx_gamma = np.power((r_x),2)*2.*h_p*nu_opt*1e9*pabs*1.e-12*(1. + n0_gamma) 

#   Define the output structure and return
    struct = {'nth':nth,'nqp':nqp,'tau_qp':tau_qp,'s1':s1,'s2':s2,'xr':xr,'Qi_inv':Qi_inv,'r_x':r_x,'r_qinv':r_qinv,'sxx_gr':sxx_gr,'sxx_gr0':sxx_gr0,'sxx_gamma':sxx_gamma}

    return struct
#***********************************************************************
##%%
## PLOTS
#t= np.arange(0.001,0.5,0.01)
#colors = ['black','m','red','coral','orange','gold','y','lightgreen','g','cyan','dodgerblue','darkorchid','hotpink']
#yLabel = ['nqp','tau_qp']#,'xr','Qi_inv','r_x','r_qinv','sxx_gr','sxx_gr0','sxx_gamma']
#power = [0,1e-15,1e-14,1e-13,1e-12,1e-11,1e-10]
#
#fig = plt.figure(1)
#fig.set_size_inches(12.5,8.5,forward=True)
#plt.clf()
#p1 = fig.add_subplot(331) # Semilog nqp vs temp
#p2 = fig.add_subplot(332) # Semilog tau_qp vs temp
#p3 = fig.add_subplot(333) 
#p4 = fig.add_subplot(334) 
#p5 = fig.add_subplot(335)
#p6 = fig.add_subplot(336)
#
#powerlevel = 'Optical Power ' + str(power[1])
#        
#p1.semilogy(t,kid_responsivity(temp=t,pabs=power[1])['nqp'],color=colors[1],label=powerlevel)
#p1.set_xlabel('Temp [K]')
#p1.set_ylabel('QP Density [QP per um^3]')
#        
#p2.semilogy(t,kid_responsivity(temp=t,pabs=power[1])['tau_qp'],color=colors[1])
#p2.set_xlabel('Temp [K]')
#p2.set_ylabel(r'$\tau_{qp}$ [s]')
#        
#p3.plot(t,kid_responsivity(temp=t,pabs=power[1])['Qi_inv'],color=colors[1])
#p3.set_xlabel('Temp [K]')
#p3.set_ylabel(r'$1/Q_i$')
#
#p4.plot(t,kid_responsivity(temp=t,pabs=power[1])['r_x'],color=colors[1])
#p4.set_xlabel('Temp [K]')
#p4.set_ylabel(r'$\Delta f/f_0$ R_x')
#
#
#p5.plot(t,kid_responsivity(temp=t,pabs=power[1])['xr'],color=colors[1])
#p5.set_xlabel('Temp [K]')
#p5.set_ylabel(r'$\Delta f/f_0$')
#
#
##p6.semilogy(t,kid_responsivity(temp=t,pabs=power[j])['sxx_gr'],color=colors[j])
##p6.set_xlabel('Temp [K]')
#        
#handles, labels = p1.get_legend_handles_labels()
#fig.legend(handles,labels,bbox_to_anchor=(.95,.6),borderaxespad=0.)
#fig.tight_layout()
#fig.show()










 