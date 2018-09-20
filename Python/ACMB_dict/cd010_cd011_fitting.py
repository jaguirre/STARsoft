# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 23:10:29 2018

@author: Alyssa
"""

import KID_model_functions_unitless as kids
import fitting_KID_model_functions as fitkids
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy import units as u


#test case:
alpha = 0.73*u.dimensionless_unscaled
f = 245*u.MHz
Tstage = 0.215*u.K
Tstage = np.linspace(0.2,0.4,15)*u.K
Tc = 1.39*u.K
TBB = 6.0*u.K
V = 76*np.power(u.micron,3)
n_star = 1318*(np.power(u.micron,-3))
tau_max = 35*u.microsecond
eta_pb = 0.57
nu_opt = (350*u.micron).to(u.GHz,equivalencies=u.spectral())
trans=1
eta_opt = 0.17*u.dimensionless_unscaled
N0=1.72e10*np.power(u.micron,-3)*np.power(u.eV,-1)

#%%
#
#TstageLTD,xLTD,QinvLTD = np.loadtxt('LTDData_x_invQi_vs_Tstage.txt',unpack=True,comments=';')
#TstageLTD*=u.K

dark2_Tstage,dark2_xavg,dark2_xerr,dark2_Qravg,dark2_Qrerr,dark2_Sxxavg,dark2_Sxxerr = np.loadtxt('cd011_res2_corr.csv',delimiter='\t',skiprows=1,unpack=True)
dark2_Tstage*=u.K
dark2_f = 336.91*u.MHz

xtest = kids.xMB(alpha,dark2_f,dark2_Tstage,Tc,5.9*u.K,V,n_star,tau_max,eta_pb,eta_opt=1,trans=0,N0=N0)
Qinvtest = kids.QinvMB(alpha,dark2_f,dark2_Tstage,Tc,5.9*u.K,V,n_star,tau_max,eta_pb,eta_opt=1,trans=0,N0=N0)
data = [dark2_Tstage,dark2_f]#[(dark2_Tstage.to(u.K)).value,(dark2_f.to(u.Hz)).value]
xsig = dark2_xerr
Qsig = dark2_Qrerr

data = [[0.215,0.25,0.275,0.3,0.325,0.35]*u.K,336.91*u.MHz]
xfitopt,xfitcov = curve_fit(fitkids.x_dark_fit,data,dark2_xavg,p0=(.73,1.39,(xtest[0]/dark2_f).value),sigma=xsig,bounds=([0,.5,(min(dark2_xavg)/f).value],[1,2.5,0]))
Qfitopt,Qfitcov = curve_fit(fitkids.Qinv_dark_fit,data,dark2_Qravg,p0=(.73,1.39,Qinvtest[0].value),bounds=([0,.5,0],[1,2.5,max(dark2_Qravg)]),sigma=Qsig)
#
#cdata = np.concatenate((dark2_xavg,dark2_Qravg))
#csigma = np.concatenate((xsig,Qsig))
#p0=(.73,1.39,(xtest[0]/f).value,Qinvtest[0].value)
#cfitopt,cfitcov = curve_fit(fitkids.x_Qinv_dark_simulfit,data,cdata,sigma=csigma,p0=p0,bounds=([0,.5,(min(xLTD)/f).value,0],[1,2.5,(-min(xLTD)/f).value,max(QinvLTD)]))
#
#plt.close('all')
#f3 = plt.figure('Figure 3')
#p1 = f3.add_subplot(121)
#p1.plot(dark2_Tstage,dark2_xavg,'ks')
#p1.plot(dark2_Tstage,xtest,'y.')
#p1.plot(TstageLTD,fitkids.x_dark_fit(data,*xfitopt),'cx')
#p1.plot(TstageLTD,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[0:len(TstageLTD)],'r-')
#
#p2 = f3.add_subplot(122)
#p2.plot(TstageLTD,QinvLTD,'ks')
#p2.plot(TstageLTD,Qinvtest,'g.')
#p2.plot(TstageLTD,fitkids.Qinv_dark_fit(data,*Qfitopt),'cx')
#p2.plot(TstageLTD,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[len(TstageLTD):],'r-')
#
##%%
#
##Pinc2LTD,x2LTD = np.loadtxt('LTDData_x_vs_Pinc.txt',unpack=True,comments=';')
##Pinc2LTD*=u.pW
##x2LTD*=u.dimensionless_unscaled
##
##PincLTD,SxxLTD = np.loadtxt('LTDData_Sxx_vs_Pinc.txt',unpack=True,comments=';')
##PincLTD*=u.pW
##SxxLTD*=np.power(u.Hz,-1)
##
#opt1_TBB,opt1_xavg,opt1_xerr,opt1_Qravg,opt1_Qrerr,opt1_Sxxavg,opt1_Sxxerr = np.loadtxt('cd010_res1_corr.csv',delimiter=',',skiprows=1,unpack=True)
#
#opt_Pinc = kids.TBB_to_Pinc(opt1_TBB)
#Tstage0 = .215*u.K
#f1 = 324.56*u.MHz
#
#
#x2test = kids.xMB(alpha,f1,Tstage0,Tc,opt1_TBB,V,n_star,tau_max,eta_pb,trans=1,eta_opt=.17,N0=N0)
#sxx2test = kids.Sxx(alpha,f1,Tstage0,Tc,opt1_TBB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt=0.17,trans=1,N0=N0)
#
#data2 = [opt1_TBB,alpha,f1,Tstage0,Tc,V,eta_pb*u.dimensionless_unscaled,nu_opt,trans,N0]
#p0x = (n_star.value,tau_max.value,eta_opt,(x2test[1]/f).value)
#boundsx = ([0,0,0,(min(x2_interp)/f).value],[np.inf,1e5,1,-(min(x2_interp)/f).value])
#x2sig = 0.05*x2_interp[1]*np.ones_like(x2_interp)
#
#p0Sxx = (n_star.value,tau_max.value,eta_opt,(Sxx_interp[0]/2).value)
#boundsSxx = ([0,0,0,0],[np.inf,1e5,1,(Sxx_interp[0]).value])
#Sxxsig = 0.05*Sxx_interp
#
#x2fitopt,x2fitcov = curve_fit(fitkids.x_opt_fit,data2,x2_interp,sigma=x2sig,p0=p0x,bounds=boundsx)
#Sxxfitopt,Sxxfitcov = curve_fit(fitkids.Sxx_fit,data2,Sxx_interp,sigma=Sxxsig,p0=p0Sxx,bounds=boundsSxx)
#
#opt_data = np.concatenate((x2_interp,Sxx_interp))
#opt_sigma = np.concatenate((x2sig,Sxxsig))
#opt_p0 =(n_star.value,tau_max.value,eta_opt,(x2test[1]/f).value,(Sxx_interp[0]/2).value)
#opt_bounds = ([0,0,0,(min(x2_interp)/f).value,0],[np.inf,1e5,1,-(min(x2_interp)/f).value,(Sxx_interp[0]).value])
#opt_fitopt,opt_fitcov = curve_fit(fitkids.x_Sxx_opt_simulfit,data2,opt_data,sigma=opt_sigma,p0=opt_p0,bounds=opt_bounds)
#
#
##plt.close('all')
#f5 = plt.figure('Figure 5')
#p3 = f5.add_subplot(121)
#p3.plot(opt_Pinc,opt1_xavg,'ks')
#p3.plot(Pinc,x2_interp,'ks')
#p3.plot(opt_Pinc,x2test,'g.')
#p3.plot(Pinc,fitkids.x_opt_fit(data2,*x2fitopt),'cx')
#p3.plot(Pinc,fitkids.x_Sxx_opt_simulfit(data2,*opt_fitopt)[0:len(Pinc)],'r-')
#
#p4 = f5.add_subplot(122)
#p4.plot(PincLTD,SxxLTD,'ys')
#p4.plot(Pinc,Sxx_interp,'ks')
#p4.plot(Pinc,sxx2test,'g.')
#p4.plot(Pinc,fitkids.Sxx_fit(data2,*Sxxfitopt),'cx')
#p4.plot(Pinc,fitkids.x_Sxx_opt_simulfit(data2,*opt_fitopt)[len(Pinc):],'r-')
#
#
#
#
