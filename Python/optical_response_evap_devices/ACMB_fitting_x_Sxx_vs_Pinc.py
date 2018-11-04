# -*- coding: utf-8 -*-
"""
Created on Sun Sep 16 16:46:46 2018

@author: Alyssa
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 23:10:29 2018
@author: Alyssa
"""

import KID_model_functions as kids
import fitting_KID_model_functions as fitkids
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy import units as u


#test case:
alpha = 0.73*u.dimensionless_unscaled
f = 330*u.MHz
Tstage = 0.215*u.K
#Tstage = np.linspace(0.2,0.4,15)*u.K
Tc = 1.39*u.K
#TBB = 6.0*u.K
V = 57*np.power(u.micron,3)
n_star = 1318*(np.power(u.micron,-3))
tau_max = 35*u.microsecond
eta_pb = 0.57
nu_opt = (350*u.micron).to(u.GHz,equivalencies=u.spectral())
trans=1
eta_opt = 0.17*u.dimensionless_unscaled
N0=1.72e10*np.power(u.micron,-3)*np.power(u.eV,-1)

#%%
#TstageLTD,xLTD,QinvLTD = np.loadtxt('LTDData_x_invQi_vs_Tstage.txt',unpack=True,comments=';')
#TstageLTD*=u.K
#
##dark2_Tstage,dark2_xavg,dark2_xerr,dark2_Qravg,dark2_Qrerr,dark2_Sxxavg,dark2_Sxxerr = np.loadtxt('cd011_res2_corr.txt',unpack=True,skiprows=1)
##TstageLTD = dark2_Tstage
##TstageLTD*=u.K
#
##xLTD = dark2_xavg
##QinvLTD = dark2_Qravg
#
#xtest = kids.xMB(alpha,f,TstageLTD,Tc,6*u.K,V,n_star,tau_max,eta_pb,eta_opt=1,trans=0,N0=N0)
#Qinvtest = kids.QinvMB(alpha,f,TstageLTD,Tc,6*u.K,V,n_star,tau_max,eta_pb,eta_opt=1,trans=0,N0=N0)
#data = [TstageLTD,f]
#xsig = 0.05*xLTD[1]*np.ones_like(xLTD)
#Qsig = 0.05*QinvLTD
#
#xfitopt,xfitcov = curve_fit(fitkids.x_dark_fit,data,xLTD,p0=(.73,1.39,(xtest[0]/f).value),sigma=xsig,bounds=([0,.5,(min(xLTD)/f).value],[1,2.5,0]))
#Qfitopt,Qfitcov = curve_fit(fitkids.Qinv_dark_fit,data,QinvLTD,p0=(.73,1.39,Qinvtest[0]),sigma=Qsig,bounds=([0,.5,0],[1,2.5,max(QinvLTD)]))
#
#cdata = np.concatenate((xLTD,QinvLTD))
#csigma = np.concatenate((xsig,Qsig))
#p0=(.73,1.39,(xtest[0]/f).value,Qinvtest[0].value)
#cfitopt,cfitcov = curve_fit(fitkids.x_Qinv_dark_simulfit,data,cdata,sigma=csigma,p0=p0,bounds=([0,.5,(min(xLTD)/f).value,0],[1,2.5,(-min(xLTD)/f).value,max(QinvLTD)]))
#
#plt.close('all')
#f3 = plt.figure('Figure 3')
#p1 = f3.add_subplot(121)
#p1.plot(TstageLTD,xLTD,'ks')
#p1.plot(TstageLTD,xtest,'g.')
#p1.plot(TstageLTD,fitkids.x_dark_fit(data,*xfitopt),'cx')
#p1.plot(TstageLTD,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[0:len(TstageLTD)],'r-')
#
#p2 = f3.add_subplot(122)
#p2.plot(TstageLTD,QinvLTD,'ks')
#p2.plot(TstageLTD,Qinvtest,'g.')
#p2.plot(TstageLTD,fitkids.Qinv_dark_fit(data,*Qfitopt),'cx')
#p2.plot(TstageLTD,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[len(TstageLTD):],'r-')

#%%

Pinc2LTDi,x2LTD,x2sig = np.loadtxt('comb_evap_x_vs_Pinc.txt',unpack=True)
Pinc2LTDi*=u.pW
x2LTD*=u.dimensionless_unscaled
indx2 = np.argsort(Pinc2LTDi)
indx22 = np.arange(0,8)
Pinc2LTD = (Pinc2LTDi[indx2])[indx22]
x2LTD = (x2LTD[indx2])[indx22]
x2sig = (x2sig[indx2])[indx22]

PincLTDi,SxxLTD,Sxxsig = np.loadtxt('comb_evap_Sxx_vs_Pinc.txt',unpack=True)
PincLTDi*=u.pW
SxxLTD*=np.power(u.Hz,-1)
indx = np.argsort(PincLTDi)
indx1 = np.arange(0,9)
PincLTD = (PincLTDi[indx])[indx1]
SxxLTD = (SxxLTD[indx])[indx1]
Sxxsig = (Sxxsig[indx])[indx1]

Tstage0 = .215*u.K
f = 330*u.MHz

# reverse engineering the BB temps
temptemps = np.linspace(0,20,1000)*u.K
temppincs = kids.TBB_to_Pinc(temptemps)
xTBB = [temptemps[np.max(np.where(temppincs<pn))].value for pn in Pinc2LTD]*u.K # goes with x data
sxxTBB = np.zeros(len(PincLTD))*u.K # goes with Sxx data
for ix,pn in enumerate(PincLTD):
    if pn!=0:
        sxxTBB[ix] = temptemps[np.max(np.where(temppincs<pn))]

    
#T_BB = np.linspace(6,7.75,15)*u.K
#Pinc = kids.TBB_to_Pinc(T_BB)
#x2_interp = np.interp(Pinc,Pinc2LTD,x2LTD)
#Sxx_interp = np.interp(Pinc,PincLTD,SxxLTD)*np.power(u.Hz,-1)
x2_interp = x2LTD
Sxx_interp = SxxLTD

x2test = kids.xMB(alpha,f,Tstage0,Tc,xTBB,V,n_star,tau_max,eta_pb,trans=1,eta_opt=.17,N0=N0)
sxx2test = kids.Sxx(alpha,f,Tstage0,Tc,sxxTBB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt=0.17,trans=1,N0=N0)

xdata = [xTBB,alpha,f,Tstage0,Tc,V,eta_pb*u.dimensionless_unscaled,nu_opt,trans,N0]
p0x = (n_star.value,tau_max.value,eta_opt,0)
boundsx = ([0,0,0,-1],[np.inf,1e5,1,1])
#boundsx = ([0,0,0,(min(x2_interp)/f).value],[np.inf,1e5,1,-(min(x2_interp)/f).value])
#x2sig = 0.05*x2_interp[1]*np.ones_like(x2_interp)

sxxdata = [sxxTBB,alpha,f,Tstage0,Tc,V,eta_pb*u.dimensionless_unscaled,nu_opt,trans,N0]
p0Sxx = (n_star.value,tau_max.value,eta_opt,(Sxx_interp[0]/2).value)
boundsSxx = ([0,0,0,0],[np.inf,1e5,1,(Sxx_interp[0]).value])
#Sxxsig = 0.05*Sxx_interp

x2fitopt,x2fitcov = curve_fit(fitkids.x_opt_fit,xdata,x2_interp,sigma=x2sig,p0=p0x,bounds=boundsx)
Sxxfitopt,Sxxfitcov = curve_fit(fitkids.Sxx_fit,sxxdata,Sxx_interp,sigma=Sxxsig,p0=p0Sxx,bounds=boundsSxx)
#
opt_data = np.concatenate((x2_interp,Sxx_interp))
opt_sigma = np.concatenate((x2sig,Sxxsig))
#opt_p0 =(n_star.value,tau_max.value,eta_opt,(x2test[1]/f).value,(Sxx_interp[0]/2).value)
opt_p0 =(Sxxfitopt[0],Sxxfitopt[1],eta_opt,(x2test[1]/f).value,(Sxx_interp[0]/2).value)
opt_bounds = ([0,0,0,-1,0],[np.inf,1e5,1,1,(Sxx_interp[0]).value])
#opt_bounds = ([0,0,0,(min(x2_interp)/f).value,0],[np.inf,1e5,1,-(min(x2_interp)/f).value,(Sxx_interp[0]).value])
opt_fitopt,opt_fitcov = curve_fit(fitkids.x_Sxx_opt_simulfit,[xdata,sxxdata],opt_data,sigma=opt_sigma,p0=opt_p0,bounds=opt_bounds)

#%%
plt.close('all')
f5 = plt.figure('Figure 5')
p3 = f5.add_subplot(121)
p3.errorbar(Pinc2LTD.value,x2LTD,yerr=x2sig,fmt='s',markerfacecolor='w',markeredgecolor='k',ecolor='k')
#p3.plot(kids.TBB_to_Pinc(xTBB).value,x2_interp,'ko',markersize=10)
#p3.plot(Pinc,x2test,'g.')
p3.plot(kids.TBB_to_Pinc(xTBB),fitkids.x_opt_fit(xdata,*x2fitopt),':',color='dodgerblue',linewidth=3.0)
p3.plot(kids.TBB_to_Pinc(xTBB),fitkids.x_Sxx_opt_simulfit([xdata,sxxdata],*opt_fitopt)[0:len(xTBB)],'r-',linewidth=2.0)
#p3.plot(kids.TBB_to_Pinc(xTBB),fitkids.x_opt_fit(xdata,Sxxfitopt[0],Sxxfitopt[1],.6,x2fitopt[-1]),'g:')
p3.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
p3.set_xlabel(r'$P_{inc}$ (pW)')
p3.set_ylabel('x = df/f')

p4 = f5.add_subplot(122,sharex=p3)
p4.errorbar(PincLTD.value,SxxLTD.value,yerr=Sxxsig,fmt='s',markerfacecolor='w',markeredgecolor='k',ecolor='k',label='data')
#p4.plot(kids.TBB_to_Pinc(sxxTBB),Sxx_interp,'ko',markersize=10)
#p4.plot(Pinc,sxx2test,'g.')
p4.plot(kids.TBB_to_Pinc(sxxTBB),fitkids.Sxx_fit(sxxdata,*Sxxfitopt),':',color='dodgerblue',linewidth=3.0,label='indiv. fits')
p4.plot(kids.TBB_to_Pinc(sxxTBB),fitkids.x_Sxx_opt_simulfit([xdata,sxxdata],*opt_fitopt)[len(xTBB):],'r-',linewidth=2.0,label='simul. fits')
#p4.plot(kids.TBB_to_Pinc(sxxTBB),fitkids.Sxx_fit(sxxdata,Sxxfitopt[0],Sxxfitopt[1],.6,Sxxfitopt[-1]),'g:')
p4.set_xlabel(r'$P_{inc}$ (pW)')
p4.set_ylabel(r'$S_{xx}$ (Hz$^{-1}$)')
p4.legend(loc='lower right')

f5.suptitle('Evap devices combined optical data')