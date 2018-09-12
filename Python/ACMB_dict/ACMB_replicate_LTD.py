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
TstageLTD,xLTD,QinvLTD = np.loadtxt('LTDData_x_invQi_vs_Tstage.txt',unpack=True,comments=';')
TstageLTD*=u.K

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
Sxx_interp = np.interp(Pinc,PincLTD,SxxLTD)

x2test = kids.xMB(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,trans=1,eta_opt=.17,N0=N0)
sxx2test = kids.Sxx(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt=0.17,trans=1,N0=N0)

data2 = [alpha,f,Tstage0,Tc,T_BB,V,eta_pb*u.dimensionless_unscaled,nu_opt,trans*u.dimensionless_unscaled,N0]
data3 = [np.asarray_chkfinite(da) for da in data2]
    
data2 = [f,alpha]
print(np.asarray_chkfinite(data2))
x2sig = 0.05*x2_interp[1]*np.ones_like(x2_interp)
Sxxsig = 0.05*Sxx_interp

#xfitopt,xfitcov = curve_fit(fitkids.x_dark_fit,data,xLTD,p0=(.73,1.39,(xtest[0]/f).value),sigma=xsig,bounds=([0,.5,(min(xLTD)/f).value],[1,2.5,0]))
x2fitopt,x2fitcov = curve_fit(fitkids.x_opt_fit,data3,x2_interp)#,sigma=x2sig)

#plt.close('all')
f5 = plt.figure('Figure 5')
p3 = f5.add_subplot(121)
p3.plot(Pinc2LTD,x2LTD,'ks')
p3.plot(kids.TBB_to_Pinc(T_BB),x2test,'g.')

p4 = f5.add_subplot(122)
p4.plot(PincLTD,SxxLTD,'ks')
p4.plot(kids.TBB_to_Pinc(T_BB),sxx2test,'g.')




