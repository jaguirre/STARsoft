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
#TstageLTD,xLTD,QinvLTD = np.loadtxt('LTDData_x_invQi_vs_Tstage.txt',unpack=True,comments=';')
#TstageLTD*=u.K

dark2_Tstage,dark2_xavg,dark2_xerr,dark2_Qravg,dark2_Qrerr,dark2_Sxxavg,dark2_Sxxerr = np.loadtxt('cd011_res2_corr.txt',unpack=True,skiprows=1)
TstageLTD = dark2_Tstage
TstageLTD*=u.K

xLTD = dark2_xavg
QinvLTD = dark2_Qravg

xtest = kids.xMB(alpha,f,TstageLTD,Tc,6*u.K,V,n_star,tau_max,eta_pb,eta_opt=1,trans=0,N0=N0)
Qinvtest = kids.QinvMB(alpha,f,TstageLTD,Tc,6*u.K,V,n_star,tau_max,eta_pb,eta_opt=1,trans=0,N0=N0)
data = [TstageLTD,f]
xsig = 0.05*xLTD[1]*np.ones_like(xLTD)
Qsig = 0.05*QinvLTD

xfitopt,xfitcov = curve_fit(fitkids.x_dark_fit,data,xLTD,p0=(.73,1.39,(xtest[0]/f).value),sigma=xsig,bounds=([0,.5,(min(xLTD)/f).value],[1,2.5,0]))
Qfitopt,Qfitcov = curve_fit(fitkids.Qinv_dark_fit,data,QinvLTD,p0=(.73,1.39,Qinvtest[0]),sigma=Qsig,bounds=([0,.5,0],[1,2.5,max(QinvLTD)]))

cdata = np.concatenate((xLTD,QinvLTD))
csigma = np.concatenate((xsig,Qsig))
p0=(.73,1.39,(xtest[0]/f).value,Qinvtest[0].value)
cfitopt,cfitcov = curve_fit(fitkids.x_Qinv_dark_simulfit,data,cdata,sigma=csigma,p0=p0,bounds=([0,.5,(min(xLTD)/f).value,0],[1,2.5,(-min(xLTD)/f).value,max(QinvLTD)]))

plt.close('all')
f3 = plt.figure('Figure 3')
p1 = f3.add_subplot(121)
p1.plot(TstageLTD,xLTD,'ks')
p1.plot(TstageLTD,xtest,'g.')
p1.plot(TstageLTD,fitkids.x_dark_fit(data,*xfitopt),'cx')
p1.plot(TstageLTD,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[0:len(TstageLTD)],'r-')

p2 = f3.add_subplot(122)
p2.plot(TstageLTD,QinvLTD,'ks')
p2.plot(TstageLTD,Qinvtest,'g.')
p2.plot(TstageLTD,fitkids.Qinv_dark_fit(data,*Qfitopt),'cx')
p2.plot(TstageLTD,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[len(TstageLTD):],'r-')

#%%

Pinc2LTD,x2LTD = np.loadtxt('LTDData_x_vs_Pinc.txt',unpack=True,comments=';')
Pinc2LTD*=u.pW
x2LTD*=u.dimensionless_unscaled

PincLTD,SxxLTD = np.loadtxt('LTDData_Sxx_vs_Pinc.txt',unpack=True,comments=';')
PincLTD*=u.pW
SxxLTD*=np.power(u.Hz,-1)

Tstage0 = .215*u.K
f = 245*u.MHz

# reverse engineering the BB temps
T_BB = np.linspace(6,11,15)*u.K
Pinc = kids.TBB_to_Pinc(T_BB)
x2_interp = np.interp(Pinc,Pinc2LTD,x2LTD)
Sxx_interp = np.interp(Pinc,PincLTD,SxxLTD)*np.power(u.Hz,-1)

x2test = kids.xMB(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,trans=1,eta_opt=.17,N0=N0)
sxx2test = kids.Sxx(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt=0.17,trans=1,N0=N0)

data2 = [T_BB,alpha,f,Tstage0,Tc,V,eta_pb*u.dimensionless_unscaled,nu_opt,trans,N0]
p0x = (n_star.value,tau_max.value,eta_opt,(x2test[1]/f).value)
boundsx = ([0,0,0,(min(x2_interp)/f).value],[np.inf,1e5,1,-(min(x2_interp)/f).value])
x2sig = 0.05*x2_interp[1]*np.ones_like(x2_interp)

p0Sxx = (n_star.value,tau_max.value,eta_opt,(Sxx_interp[0]/2).value)
boundsSxx = ([0,0,0,0],[np.inf,1e5,1,(Sxx_interp[0]).value])
Sxxsig = 0.05*Sxx_interp

x2fitopt,x2fitcov = curve_fit(fitkids.x_opt_fit,data2,x2_interp,sigma=x2sig,p0=p0x,bounds=boundsx)
Sxxfitopt,Sxxfitcov = curve_fit(fitkids.Sxx_fit,data2,Sxx_interp,sigma=Sxxsig,p0=p0Sxx,bounds=boundsSxx)

opt_data = np.concatenate((x2_interp,Sxx_interp))
opt_sigma = np.concatenate((x2sig,Sxxsig))
opt_p0 =(n_star.value,tau_max.value,eta_opt,(x2test[1]/f).value,(Sxx_interp[0]/2).value)
opt_bounds = ([0,0,0,(min(x2_interp)/f).value,0],[np.inf,1e5,1,-(min(x2_interp)/f).value,(Sxx_interp[0]).value])
opt_fitopt,opt_fitcov = curve_fit(fitkids.x_Sxx_opt_simulfit,data2,opt_data,sigma=opt_sigma,p0=opt_p0,bounds=opt_bounds)


#plt.close('all')
f5 = plt.figure('Figure 5')
p3 = f5.add_subplot(121)
p3.plot(Pinc2LTD,x2LTD,'ys')
p3.plot(Pinc,x2_interp,'ks')
p3.plot(Pinc,x2test,'g.')
p3.plot(Pinc,fitkids.x_opt_fit(data2,*x2fitopt),'cx')
p3.plot(Pinc,fitkids.x_Sxx_opt_simulfit(data2,*opt_fitopt)[0:len(Pinc)],'r-')

p4 = f5.add_subplot(122)
p4.plot(PincLTD,SxxLTD,'ys')
p4.plot(Pinc,Sxx_interp,'ks')
p4.plot(Pinc,sxx2test,'g.')
p4.plot(Pinc,fitkids.Sxx_fit(data2,*Sxxfitopt),'cx')
p4.plot(Pinc,fitkids.x_Sxx_opt_simulfit(data2,*opt_fitopt)[len(Pinc):],'r-')