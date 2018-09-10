# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 23:10:29 2018

@author: Alyssa
"""

import KID_model_functions as kids
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy import units as u


#test case:
alpha = 0.8*u.dimensionless_unscaled
f = 250*u.MHz
Tstage = 0.215*u.K
Tstage = np.linspace(0.2,0.4,15)*u.K
Tc = 1.4*u.K
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

xtest = kids.xMB(alpha,f,TstageLTD,Tc,6*u.K,V,n_star,tau_max,eta_pb,0,eta_opt=1,N0=N0)
Qinvtest = kids.QinvMB(alpha,f,TstageLTD,Tc,0,V,n_star,tau_max,eta_pb,0,N0)
data = [TstageLTD,f]
xsig = 0.05*xLTD[1]*np.ones_like(xLTD)
Qsig = 0.05*QinvLTD

a = kids.xMB(alpha,f,TstageLTD,Tc,TBB=1,V=1,n_star=1,tau_max=1,eta_pb=1,trans=0,eta_opt=1,N0=N0)
xfitopt,xfitcov = curve_fit(kids.x_dark_fit,data,xLTD,p0=(.73,1.39,xtest[0]),sigma=xsig,bounds=([0,.5,min(xLTD)],[1,2.5,0]))
Qfitopt,Qfitcov = curve_fit(kids.Qinv_dark_fit,data,QinvLTD,p0=(.73,1.39,Qinvtest[0]),sigma=Qsig,bounds=([0,.5,0],[1,2.5,max(QinvLTD)]))

cdata = np.concatenate((xLTD,QinvLTD))
csigma = np.concatenate((xsig,Qsig))
p0=(.73,1.39,xtest[0],Qinvtest[0])
cfitopt,cfitcov = curve_fit(kids.x_Qinv_dark_simulfit,data,cdata,sigma=csigma,p0=p0,bounds=([0,.5,min(xLTD),0],[1,2.5,-min(xLTD),max(QinvLTD)]))

plt.close('all')
f3 = plt.figure('Figure 3')
p1 = f3.add_subplot(121)
p1.plot(TstageLTD,xLTD,'ks')
p1.plot(TstageLTD,xtest,'g.')
p1.plot(TstageLTD,kids.x_dark_fit(data,*xfitopt),'cx')
p1.plot(TstageLTD,kids.x_Qinv_dark_simulfit(data,*cfitopt)[0:len(TstageLTD)],'r-')

p2 = f3.add_subplot(122)
p2.plot(TstageLTD,QinvLTD,'ks')
p2.plot(TstageLTD,Qinvtest,'g.')
p2.plot(TstageLTD,kids.Qinv_dark_fit(data,*Qfitopt),'cx')
p2.plot(TstageLTD,kids.x_Qinv_dark_simulfit(data,*cfitopt)[len(TstageLTD):],'r-')

#%%

Pinc2LTD,x2LTD = np.loadtxt('LTDData_x_vs_Pinc.txt',unpack=True,comments=';')
Pinc2LTD*=u.pW
x2LTD*=u.dimensionless_unscaled

PincLTD,SxxLTD = np.loadtxt('LTDData_Sxx_vs_Pinc.txt',unpack=True,comments=';')
PincLTD*=u.pW
SxxLTD*=np.power(u.Hz,-1)

Tstage0 = .215*u.K
T_BB = np.linspace(6,12,15)*u.K
x2test = kids.xMB(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,trans=1,eta_opt=.17,N0=N0)

sxx2test = kids.Sxx(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,trans=1,eta_opt=1,N0=N0)
data = [alpha,f,Tstage,TBB,V,eta_pb,nu_opt,trans,N0]

plt.close('all')
f5 = plt.figure('Figure 5')
p3 = f5.add_subplot(121)
p3.plot(Pinc2LTD,x2LTD,'ks')
p3.plot(kids.TBB_to_Pinc(T_BB),x2test,'g.')

p4 = f5.add_subplot(122)
p4.plot(PincLTD,SxxLTD,'ks')
p4.plot(kids.TBB_to_Pinc(T_BB),sxx2test,'g.')


