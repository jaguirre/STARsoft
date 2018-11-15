#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 10:45:05 2018

@author: starfire
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
#import KID_model_functions as kids
#import fitting_KID_model_functions as fitkids

res = 0
plt.close('all')
fig,(p1,p2) = plt.subplots(2,1,sharex=True,num=('Res '+str(res)))

cd010 = {'CD':'CD010',
         'trans': 1,
         'color': 'blue'}

cd012 = {'CD':'CD012',
         'trans': 0.03,
         'color': 'red'}

cd013 = {'CD':'CD013',
         'trans': 1,
         'color': 'purple'}

for cool in [cd010,cd012,cd013]:
    datafile = cool['CD']+'_reduced_Res'+str(res)+'.csv'
    
    TBBlist,f0list,xlist,Sxxlist = np.loadtxt(datafile,delimiter=',',unpack=True)
    cool['TBBlist'] = TBBlist
    cool['f0list'] = f0list
    cool['xlist'] = xlist
    cool['Sxxlist'] = Sxxlist
    
    TBBavg = list(set(TBBlist))
    f0avg = []
    f0err = []
    xavg = []
    xerr = []
    Sxxavg = []
    Sxxerr = []
    
    for temp in TBBavg:
        inds = np.where(TBBlist==temp)
        f0avg.append(np.average(f0list[inds]))        
        f0err.append(np.std(f0list[inds])/np.sqrt(len(inds)))
        xavg.append(np.average(xlist[inds]))
        xerr.append(np.std(xlist[inds])/np.sqrt(len(inds)))
        Sxxavg.append(np.average(Sxxlist[inds]))
        Sxxerr.append(np.std(Sxxlist[inds])/np.sqrt(len(inds)))
    
    cool['TBBavg'] = TBBavg
    cool['f0avg'] = f0avg
    cool['xavg'] = xavg
    cool['Sxxavg'] = Sxxavg
    
    cool['f0err'] = f0err
    cool['xerr'] = xerr
    cool['Sxxerr'] = Sxxerr


cd011 = {'CD':'CD011',
         'trans':0,
         'color':'green'}
datafile = cd011['CD']+'_reduced_Res'+str(res)+'.csv'
Tstagelist,f0list,xlist,Sxxlist = np.loadtxt(datafile,delimiter=',',unpack=True)
indxs = np.where(Tstagelist==.215)
f0list = f0list[indxs]
xlist = xlist[indxs]
Sxxlist = Sxxlist[indxs]
Popt = np.zeros_like(xlist)

cd011['TBBlist'] = np.zeros_like(f0list)
cd011['f0list'] = f0list
cd011['xlist'] = xlist
cd011['Sxxlist'] = Sxxlist

f0avg = []
f0err = []
xavg = []
xerr = []
Sxxavg = []
Sxxerr = []

f0avg.append(np.average(f0list))        
f0err.append(np.std(f0list)/np.sqrt(len(f0list)))
xavg.append(np.average(xlist))
xerr.append(np.std(xlist)/np.sqrt(len(xlist)))
Sxxavg.append(np.average(Sxxlist))
Sxxerr.append(np.std(Sxxlist)/np.sqrt(len(Sxxlist)))

cd011['TBBavg'] = [0]
cd011['f0avg'] = f0avg
cd011['xavg'] = xavg
cd011['Sxxavg'] = Sxxavg

cd011['f0err'] = f0err
cd011['xerr'] = xerr
cd011['Sxxerr'] = Sxxerr

#%%
''' Plot & fit! '''
pmax = 0.5*u.pW
xPtofit = []
sxxtofit = []

    
for cool in [cd010,cd011,cd012,cd013]:
    Popt = TBB_to_Pinc(cool['TBBlist'],trans=cool['trans'])
    cool['Popt'] = Popt
    Poptavg = TBB_to_Pinc(cool['TBBavg'],trans=cool['trans'])
    cool['Poptavg'] = Poptavg


    p2.plot(Popt,cool['Sxxlist'],'.',color=cool['color'],alpha=0.25)
    p2.errorbar(Poptavg.value,cool['Sxxavg'],yerr=cool['Sxxerr'],fmt='d',markersize=5.0,color=cool['color'],label=cool['CD'])
    
    
    if cool==cd011:
        sxxtofit.append(np.array((cool['TBBavg'],Poptavg.value,cool['Sxxavg'][0],cool['Sxxerr'][0])))
    else:    
        for m,n in enumerate(Poptavg): 
            if n < pmax: sxxtofit.append(np.array((cool['TBBavg'][m],Poptavg[m].value,cool['Sxxavg'][m],cool['Sxxerr'][m])))

sxxtofit = np.transpose(sxxtofit)

# fit cd012 x data to a line
z = np.polyfit(cd012['Popt'][np.where(cd012['Popt']<pmax)],cd012['xlist'][np.where(cd012['Popt']<pmax)],1)

# Use the fit and the lowest-power points from the optical cooldowns to scale the x data
x013 = cd013['Popt'][np.where(cd013['Popt']==np.min(cd013['Popt']))]
y013 = cd013['xlist'][np.where(cd013['Popt']==np.min(cd013['Popt']))]
yn013 = np.polyval(z,x013)
d013 = np.mean(y013-yn013)
yf013 = cd013['xlist']-d013

pcom = cd010['Popt'][30]
ind010 = np.where(cd010['Popt'] == pcom)
av010 = np.mean(cd010['xlist'][ind010])

ind013 = np.where(cd013['Popt'] == pcom)
av013 = np.mean(yf013[ind013])
dav = av010-av013

yf010 = cd010['xlist']-dav

cd013['xavg'] = cd013['xavg']-d013
cd013['xlist'] = yf013
cd010['xavg'] = cd010['xavg']-dav
cd010['xlist'] = yf010

for cool in [cd010,cd012,cd013]:
    
    Popt = cool['Popt']
    Poptavg = cool['Poptavg']
    
    # fit f0 vs TBB to a line
#    z = np.polyfit(Popt[np.where(Popt<pmax)],cool['xlist'][np.where(Popt<pmax)],1)
#    p1.plot(Popt,np.polyval(z,Popt)-z[1],'-',color=cool['color'])
#    p1.plot(Popt,(cool['xlist']-z[1]),'.',color=cool['color'],alpha=0.25)
    p1.errorbar(Poptavg.value,(cool['xavg']),yerr=cool['xerr'],fmt='d',markersize=5.0,color=cool['color'],label=cool['CD'])
    p1.plot(Popt,cool['xlist'],'.',color=cool['color'],alpha=0.25)
    for m,n in enumerate(Poptavg): 
        if n < pmax: xPtofit.append(np.array((cool['TBBavg'][m],Poptavg[m].value,cool['xavg'][m],cool['xerr'][m])))
    
xPtofit = np.transpose(xPtofit)


Pinc2LTD,x2LTD = np.loadtxt('../ACMB_dict/LTDData_x_vs_Pinc.txt',unpack=True,comments=';')
Pinc2LTD*=u.pW
x2LTD*=u.dimensionless_unscaled
ynLTD = np.polyval(z,Pinc2LTD[0])
dLTD = x2LTD[0]-ynLTD
x2LTD = x2LTD-dLTD

p1.plot(Pinc2LTD,x2LTD,'ko',markersize=5.0)

PincLTD,SxxLTD = np.loadtxt('../ACMB_dict/LTDData_Sxx_vs_Pinc.txt',unpack=True,comments=';')
p2.plot(PincLTD,SxxLTD,'ko',markersize=5.0,label='sputtered')


p1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
p1.set_ylabel('df/f')
#p1.legend()
p2.set_ylabel('amp sub Sxx (1/Hz)')
p2.set_xlabel(r'$P_{inc}$ (pW)')
#p2.set_xlim([-.05,0.5])
p1.legend(handles=p2.get_legend_handles_labels()[0],labels=p2.get_legend_handles_labels()[1],loc='upper right',ncol=1,fontsize=10)
#p2.legend(loc='lower right',fontsize=10)

fig.suptitle('Optical tests: sputtered vs evap')

    
#%%
fn,(pn1,pn2)=plt.subplots(2,1,sharex=True,num='fit')
#pn1.errorbar(cd012['Poptavg'].value,cd012['xavg'],yerr=cd012['xerr'],fmt='d',markersize=5.0,color='k')
pn1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pn1.set_ylabel('df/f')
pn2.set_ylabel('amp sub Sxx (1/Hz)')
pn2.set_xlabel(r'$P_{inc}$ (pW)')
pn2.set_xlim([0,.5])

pn1.errorbar(xPtofit[1],xPtofit[2],yerr=xPtofit[3],fmt='d',markersize=5.0,color='b')
pn2.errorbar(sxxtofit[1],sxxtofit[2],yerr=sxxtofit[3],fmt='d',markersize=5.0,color='b',label='evap devices')

pn1.axvline(x=0.12,linestyle=':',color='m')
pn1.axvline(x=0.2,linestyle=':',color='m')
pn2.axvline(x=0.12,linestyle=':',color='m')
pn2.axvline(x=0.2,linestyle=':',color='m')

inx = np.argsort(xPtofit[1])
xPtofit = [xPtofit[k][inx] for k in np.arange(len(xPtofit))]
inx2 = np.argsort(sxxtofit[1])
sxxtofit = [sxxtofit[k][inx2] for k in np.arange(len(sxxtofit))]

i1 = np.max(np.where(xPtofit[1]<0.1))
i2 = np.min(np.where(xPtofit[1]>0.2))
pn1.plot(xPtofit[1][i1],xPtofit[2][i1],'*',color='navy',markersize=10)
pn1.plot(xPtofit[1][i2],xPtofit[2][i2],'*',color='navy',markersize=10)

R_est = ((xPtofit[2][i2]-xPtofit[2][i1])/(u.pW*(xPtofit[1][i2]-xPtofit[1][i1]))).to(np.power(u.W,-1))

pn1.plot(Pinc2LTD[1],x2LTD[1],'ko',markersize=5)
pn1.plot(Pinc2LTD[0],x2LTD[0],'ko',markersize=5)

R_LTD = ((x2LTD[0]-x2LTD[1])/((Pinc2LTD[0]-Pinc2LTD[1]))).to(np.power(u.W,-1))

i3 = np.max(np.where(sxxtofit[1]<0.1))
i4 = np.min(np.where(sxxtofit[1]>0.2))
pn2.plot(sxxtofit[1][i3],sxxtofit[2][i3],'*',color='navy',markersize=10)
pn2.plot(sxxtofit[1][i4],sxxtofit[2][i4],'*',color='navy',markersize=10)

sxx_avg = np.power(u.Hz,-1)*(sxxtofit[2][i3]+sxxtofit[2][i4])/2
#pn2.plot((sxxtofit[1][i3]+sxxtofit[1][i4])/2,sxx_avg,'x',color='b',markersize=10)
pn2.axhline(y=sxx_avg.value,color='b',linestyle='--')

pn2.plot(PincLTD[1],SxxLTD[1],'ko',markersize=5)
pn2.plot(PincLTD[2],SxxLTD[2],'ko',markersize=5)

sxxLTD_avg = np.power(u.Hz,-1)*(SxxLTD[1]+SxxLTD[2])/2
#pn2.plot((PincLTD[1]+PincLTD[2])/2,sxxLTD_avg,'kx',markersize=8)
pn2.axhline(y=sxxLTD_avg.value,color='k',linestyle='--')

NEP_est = (np.sqrt(sxx_avg)/R_est).to(u.W*np.power(u.Hz,-.5))
NEP_LTD = (np.sqrt(sxxLTD_avg)/R_LTD).to(u.W*np.power(u.Hz,-.5))

NEP_dark = (np.sqrt(sxxtofit[2][0]*np.power(u.Hz,-1))/(xPtofit[2][i1]/(xPtofit[1][i1]*u.pW))).to(u.W*np.power(u.Hz,-.5))
NEP_LTD_dark = (np.sqrt(SxxLTD[0]*np.power(u.Hz,-1))/(x2LTD[0]/(Pinc2LTD[0]))).to(u.W*np.power(u.Hz,-.5))
#txt="Blue diamonds are combined data from the evaporated devices (cooldowns 10,11,12,13) "
#plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)
#%%
#np.savetxt('comb_evap_x_vs_Pinc.txt',np.transpose(xPtofit[1:4]),header='Popt,x=df/f,xsigma')
#np.savetxt('comb_evap_Sxx_vs_Pinc.txt',np.transpose(sxxtofit[1:4]),header='Popt,Sxx (1/Hz),Sxxsigma')

#%%
#from scipy.optimize import curve_fit
#
##test case:
#alpha = 0.73*u.dimensionless_unscaled
#f = cd012['f0list'][0]*u.Hz
#Tstage0 = 0.215*u.K
#Tc = 1.39*u.K
##TBB = 6.0*u.K
#V = 76*np.power(u.micron,3)
#n_star = 1318*(np.power(u.micron,-3))
#tau_max = 35*u.microsecond
#eta_pb = 0.57
#nu_opt = (350*u.micron).to(u.GHz,equivalencies=u.spectral())
#trans=cd012['trans']
#eta_opt = 0.17*u.dimensionless_unscaled
#N0=1.72e10*np.power(u.micron,-3)*np.power(u.eV,-1)
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
#
##Tstage0 = .215*u.K
##f = cd012['f0list'][0]*u.Hz
#
## reverse engineering the BB temps
##T_BB = np.linspace(6,11,15)*u.K
##x2_interp = np.interp(Pinc,Pinc2LTD,x2LTD)
##Sxx_interp = np.interp(Pinc,PincLTD,SxxLTD)*np.power(u.Hz,-1)
#T_BB = cd012['TBBavg']*u.K
#Pinc = TBB_to_Pinc(T_BB)
#x2_interp = cd012['xavg']
#Sxx_interp = cd012['Sxxavg']*np.power(u.Hz,-1)
#
#x2test = xMB(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,trans=1,eta_opt=.17,N0=N0)
#sxx2test = Sxx(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt=0.17,trans=1,N0=N0)
#
#data2 = [T_BB,alpha,f,Tstage0,Tc,V,eta_pb*u.dimensionless_unscaled,nu_opt,trans,N0]
#p0x = (n_star.value,tau_max.value,eta_opt,0)#(x2test[1]/f).value)
#boundsx = ([0,0,0,(min(x2_interp)/f).value],[np.inf,1e5,1,-(min(x2_interp)/f).value])
#x2sig = cd012['xerr']#0.05*x2_interp[1]*np.ones_like(x2_interp)
#
#p0Sxx = (n_star.value,tau_max.value,eta_opt,(Sxx_interp[0]/2).value)
#boundsSxx = ([0,0,0,0],[np.inf,1e5,1,(Sxx_interp[0]).value])
#Sxxsig = cd012['Sxxavg']#0.05*Sxx_interp
#
#x2fitopt,x2fitcov = curve_fit(x_opt_fit,data2,x2_interp,sigma=x2sig,p0=p0x,bounds=boundsx)
#Sxxfitopt,Sxxfitcov = curve_fit(Sxx_fit,data2,Sxx_interp,sigma=Sxxsig,p0=p0Sxx,bounds=boundsSxx)
#
#opt_data = np.concatenate((x2_interp,Sxx_interp))
#opt_sigma = np.concatenate((x2sig,Sxxsig))
#opt_p0 =(n_star.value,tau_max.value,eta_opt,0,(Sxx_interp[0]/2).value)
#opt_bounds = ([0,0,0,(min(x2_interp)/f).value,0],[np.inf,1e5,1,-(min(x2_interp)/f).value,(Sxx_interp[0]).value])
#opt_fitopt,opt_fitcov = curve_fit(x_Sxx_opt_simulfit,data2,opt_data,sigma=opt_sigma,p0=opt_p0,bounds=opt_bounds)
#
#
##plt.close('all')
#f5 = plt.figure('Figure 5')
#p3 = f5.add_subplot(121)
#p3.plot(Pinc2LTD,x2LTD,'ys')
#p3.plot(Pinc,x2_interp,'ks')
##p3.plot(Pinc,x2test,'g.')
#p3.plot(Pinc,x_opt_fit(data2,*x2fitopt),'cx')
#p3.plot(Pinc,x_Sxx_opt_simulfit(data2,*opt_fitopt)[0:len(Pinc)],'r-')
#
#p4 = f5.add_subplot(122)
#p4.plot(PincLTD,SxxLTD,'ys')
#p4.plot(Pinc,Sxx_interp,'ks')
##p4.plot(Pinc,sxx2test,'g.')
#p4.plot(Pinc,fitkids.Sxx_fit(data2,*Sxxfitopt),'cx')
#p4.plot(Pinc,fitkids.x_Sxx_opt_simulfit(data2,*opt_fitopt)[len(Pinc):],'r-')
#
#
#pn2.plot(Pinc,Sxx_interp,'ks')
##p4.plot(Pinc,sxx2test,'g.')
#pn2.plot(Pinc,Sxx_fit(data2,*Sxxfitopt),'cx')
#pn2.plot(Pinc,x_Sxx_opt_simulfit(data2,*opt_fitopt)[len(Pinc):],'r-')