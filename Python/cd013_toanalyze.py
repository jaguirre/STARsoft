# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 08:43:24 2018

@author: Alyssa
"""

from astropy import units as u

''' CD013 mixer data '''
TBB_list = u.K*[5.65,
                5.65,
                9,
                10,
                11]

Tstage_list = u.mK*[.215,
                    .215,
                    .215,
                    .215,
                    .215]

dates_list = len(TBB_list)*['20180702']

scans_list = ['noise06',
              'noise07',
              'noise10',
              'noise12',
              'noise14']

res1 = 0.33163*u.GHz
res2 = 0.34217*u.GHz
res3 = 0.35553*u.GHz
res4 = 0.36178*u.GHz
resonators = [res1,res2,res3,res4]

Qrange = 300
df = [f0/Qrange for f0 in resonators]



if (len(TBB_list)==len(Tstage_list)==len(dates_list)==len(scans_list))==False:
    print('Missed something! Input parameters are not the same length.')
#%%
#from ReadMixerData import *
#
#cool = 'CD013'
#parentdatafolder = '/scr/starfire/labdata/'
#parentoutfolder = '/scr/starfire/analysis/'
#
##for ind in [7]: #
#for ind in np.arange(0,len(scans_list)):
#    testdict = {}
#    datafolder = parentdatafolder + dates_list[ind] + '/' + scans_list[ind] + '/'
#    datescan = dates_list[ind] + '_' + scans_list[ind]
#    outfolder = parentoutfolder + cool + '/' + datescan + '/'
#    T_stage = Tstage_list[ind]
#    T_BB = TBB_list[ind]
#    #importmixerdata(testdict,T_stage,T_BB,cool,datafolder,outfolder,Pn=1,Fn=1,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5)
#    
#    importmixerfolder(testdict,T_stage,T_BB,cool,datafolder,outfolder,datescan,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5) 
#
#    # label the resonators
#    for pn in testdict.keys():
#        try:
#            for fn in testdict[pn].keys():
#                for indx,reso in enumerate(resonators):
#                    if abs(reso-testdict[pn][fn]['initial f0']*u.Hz) < df[indx]:
#                        testdict[pn][fn]['res'] = indx
#
#        except AttributeError: exit
#
#
#    savedatadict(testdict,cool,datescan,outfolder)
#    
    
#%%
import matplotlib.pyplot as plt
import numpy as np
import KID_model_functions as kids
import fitting_KID_model_functions as fitkids
from DictionaryToHDF5 import *
from dictionary_functions import *
from scipy.optimize import curve_fit
from astropy import units as u


#test case:
alpha = 0.73*u.dimensionless_unscaled
#f = 245*u.MHz
Tstage = 0.215*u.K
Tstage = np.linspace(0.2,0.4,15)*u.K
Tc = 1.39*u.K
TBB = 6.0*u.K
V = 76*np.power(u.micron,3)
n_star = 1318*(np.power(u.micron,-3))
tau_max = 35*u.microsecond
eta_pb = 0.57
nu_opt = (350*u.micron).to(u.GHz,equivalencies=u.spectral())
trans=0.03
eta_opt = 0.17*u.dimensionless_unscaled
N0=1.72e10*np.power(u.micron,-3)*np.power(u.eV,-1)

Tstage0 = 0.215*u.K

n = len(scans_list)

#for indx,res in enumerate(resonators):
indx=0
res=resonators[indx]
TBBlist = []
f0list = []
Qrinvlist = []
Sxxlist = []

for ind in np.arange(0,n):
    scan = scans_list[ind]
    testdict = hdf5_to_dict('/scr/starfire/analysis/CD013/20180702_'+scan+'/CD013_20180702_'+scan+'_datadict.hdf5')

    allpaths = dictwhere(testdict,'res',indx)    
    a = dictget(testdict,allpaths,'LB_atten')
    
    paths = [allpaths[pa] for pa in np.where(a==np.min(a))[0]]
    streampaths = [[paths[k][0],paths[k][1],'stream'] for k in np.arange(len(paths))]
    
    try:
        T_BB = dictget(testdict,paths,'T_BB')
        for T in T_BB: TBBlist.append(T)
        
        f0 = dictget(testdict,paths,'f0_fit')
        for fi in f0: f0list.append(fi)
        
        Qr = dictget(testdict,paths,'Qr_fit')
        for Q in Qr: Qrinvlist.append(1./Q)
        
        Sxx = dictget(testdict,streampaths,'amp sub Sxx white')
        for S in Sxx: Sxxlist.append(S)
    except KeyError: exit

f = f0list[0]*u.Hz    
xlist = (f0list-np.max(f0list))/np.max(f0list)
#    Sxxlist = np.array(Sxxlist)
#    x2test = kids.xMB(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,eta_opt=.17,trans=trans,N0=N0)
#    sxx2test = kids.Sxx(alpha,f,Tstage0,Tc,T_BB,V,n_star,tau_max,eta_pb,nu_opt,eta_opt=0.17,trans=trans,N0=N0)
#    
#    data2 = [TBBlist,alpha,f,Tstage0,Tc,V,eta_pb*u.dimensionless_unscaled,nu_opt,trans,N0]
#    p0x = (n_star.value,tau_max.value,eta_opt,(x2test[1]/f).value)
#    boundsx = ([0,0,0,-1],[np.inf,1e5,1,1])
##    boundsx = ([0,0,0,(min(xlist)/f).value],[np.inf,1e5,1,-(min(xlist)/f).value])
#    x2sig = 0.05*xlist[1]*np.ones_like(xlist)
#    
#    p0Sxx = (n_star.value,tau_max.value,eta_opt,(Sxxlist[0]/2))
#    boundsSxx = ([0,0,0,0],[np.inf,1e5,1,(Sxxlist[0])])
#    Sxxsig = 0.05*Sxxlist
#    
#    x2fitopt,x2fitcov = curve_fit(fitkids.x_opt_fit,data2,xlist,sigma=x2sig,p0=p0x,bounds=boundsx)
#    Sxxfitopt,Sxxfitcov = curve_fit(fitkids.Sxx_fit,data2,Sxxlist,sigma=Sxxsig,p0=p0Sxx,bounds=boundsSxx)
#    
#    opt_data = np.concatenate((x2_interp,Sxx_interp))
#    opt_sigma = np.concatenate((x2sig,Sxxsig))
#    opt_p0 =(n_star.value,tau_max.value,eta_opt,(x2test[1]/f).value,(Sxx_interp[0]/2))
#    opt_bounds = ([0,0,0,(min(x2_interp)/f).value,0],[np.inf,1e5,1,-(min(x2_interp)/f).value,(Sxx_interp[0]).value])
#    opt_fitopt,opt_fitcov = curve_fit(fitkids.x_Sxx_opt_simulfit,data2,opt_data,sigma=opt_sigma,p0=opt_p0,bounds=opt_bounds)
#    
#    fig,(p1,p2,p3) = plt.subplots(3,1,sharex=True,num=('Res '+str(indx)))
fig = plt.figure('Res ' + str(indx))
p1 = fig.add_subplot(111)
p1.plot(TBB_to_Pinc(TBBlist,trans=1),xlist,'r.')
p1.set_xlabel(r'$P_{inc}$')
p1.set_ylabel('df/f')
#    p1.plot(Tstagelist,xtest,'g--')
#    p1.plot(Tstagelist,fitkids.x_dark_fit(data,*xfitopt),'cx')
#    p1.plot(Tstagelist,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[0:len(Tstagelist)],'r-')
#    p2.plot(Tstagelist,Qrinvlist,'ko')
#    p2.plot(Tstagelist,Qinvtest,'g--')
#    p2.plot(Tstagelist,fitkids.Qinv_dark_fit(data,*Qfitopt),'cx')
#    p2.plot(Tstagelist,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[len(Tstagelist):],'r-')
#    p3.plot(Tstagelist,Sxxlist,'ko')
#    
#

np.savetxt('/scr/starfire/analysis/CD013/CD013_reduced_'+'Res'+str(indx)+'.csv',np.transpose(np.array((TBBlist,f0list,xlist,Sxxlist))),delimiter=',')

#%%
#plt.figure(100)
#colors = ['maroon','red','darkorange','green','limegreen','b','dodgerblue','m']
#
#res1 = 0.33163*u.GHz
#res2 = 0.34217*u.GHz
#res3 = 0.35553*u.GHz
#res4 = 0.36178*u.GHz
#resonators = [res1,res2,res3,res4]
#
#Qrange = 300
#df = [f0/Qrange for f0 in resonators]
#
#for res,resfreq in enumerate(resonators):
#    
#    res=0
#    resfreq=resonators[res]
#    fig,(p1,p2,p3) = plt.subplots(3,1,sharex=True,num=(str(res)+' LB'))
#    fig2,(px,psxx) = plt.subplots(2,1,sharex=True,num=(str(res)+' TBB'))
#    
#    for ind,scan in enumerate(scans_list[:-1]):
#    #ind=0
#    #scan=scans_list[ind]
#        testdict = hdf5_to_dict('/scr/starfire/analysis/CD010/20180328_'+scan+'/CD010_20180328_'+scan+'_datadict.hdf5')
#        
#        # Label each dataset by resonator
#        for pn in testdict.keys():
#            try:
#                for fn in testdict[pn].keys():
#    #                print(str(pn)+', '+str(fn))
#                    for indx,reso in enumerate(resonators):
#                        if abs(reso-testdict[pn][fn]['initial f0']*u.Hz) < df[indx]:
#                            testdict[pn][fn]['res'] = indx
#    
#            except AttributeError: exit
#        
#        ### want: P_LB vs a_nl,Qr,Qc,Qi,
#        # separate out by resonator
#    #    for res,resfreq in enumerate(resonators):
#        paths = dictwhere(testdict,'res',res)
#        streampaths = [[paths[k][0],paths[k][1],'stream'] for k in np.arange(len(paths))]
##        try:
#        p1.plot(dictget(testdict,paths,'LB_atten'),dictget(testdict,paths,'anl_fit'),'o',color=colors[ind],alpha=0.4)
#        p2.plot(dictget(testdict,paths,'LB_atten'),dictget(testdict,paths,'f0_fit'),'o',color=colors[ind],alpha=0.4)
#        p3.plot(dictget(testdict,paths,'LB_atten'),dictget(testdict,paths,'Qr_fit'),'s',color=colors[ind],alpha=0.4)
#        p3.plot(dictget(testdict,paths,'LB_atten'),dictget(testdict,paths,'Qc_fit'),'v',color=colors[ind],alpha=0.4)
#        
#        px.plot(dictget(testdict,paths,'T_BB'),dictget(testdict,paths,'f0_fit'),'ko',alpha=0.4)
#        psxx.plot(dictget(testdict,paths,'T_BB'),dictget(testdict,streampaths,'amp sub Sxx white'),'ko',alpha=0.4)
#    p1.set_xlim([25,45])
#    p1.set_ylim([-1,.25])
##    avals = list(gen_dict_extract('anl_fit',testdict))
##    pvals = list(gen_dict_extract('LB_atten',testdict))
##    plt.plot(pvals,avals,'.',color=colors[ind],alpha=0.5)
#    
##%%    
#    
#''' Resonance frequencies of resonators we studied for this device: '''
#res0 = 324.56*u.MHz
#res1 = 331.6*u.MHz
#res2 = 336.84*u.MHz
#res3 = 342.22*u.MHz
#res4 = 352.03*u.MHz
#res5 = 355.61*u.MHz
#res6 = 361.85*u.MHz
#res7 = 363.59*u.MHz
#resonators=[res0,res1,res2,res3,res4,res5,res6,res7]
#
#
#''' We'll use the range df=f0/Qrange to determine which resonator is which
#    [not the most elegant solution if there are lots of tightly-packed resonators and collisions] '''
#Qrange = 300
#df = [f0/Qrange for f0 in resonators]
#
#
#
#for cool in cdlist:
#    
#    ''' Compare each resonance frequency to the known list of resonators, 
#        use the value to index each scan by resonator number '''
#    for freq in list(gen_dict_extract('f0',cool)):
#        path = dictwhere(cool,'f0',freq)
#        for ind,res in enumerate(resonators):
#            if abs(res-freq) < df[ind]:
#                for sweep,scan in path: # hard-coding in that each freq will be indexed by cd[sweep][scan]
#                    cool[sweep][scan]['res'] = ind    
#    
#    ''' Find the maximum resonance frequency for each resonator within a cooldown,
#        set the maximum to be f00 for that resonator,cooldown '''
#    f00 = np.zeros(len(resonators))*u.MHz
#    for n,r in enumerate(resonators):
#        try: 
#            f0n = max(dictget(cool,dictwhere(cool,'res',n),'f0'))
#            f00[n] = f0n
#            
#        # if a resonator doesn't appear in this cooldown, set its f00 = -1
#        except ValueError: 
#            f00[n] = -1*u.MHz
#    cool['f00'] = f00
#    
#    ''' Use the f00 values for each resonator,cooldown to calculate x = (f-f00)/f
#        for each scan in the cooldown '''
#    for freq in list(gen_dict_extract('f0',cool)):
#        path2 = dictwhere(cool,'f0',freq)
#        reso = dictget(cool,path2,'res')
#        x = ((freq-f00[reso])/f00[reso]).to(u.dimensionless_unscaled)
#        for sweep2,scan2 in path2: # hard-coding in that each freq will be indexed by cd[sweep][scan]
#            cool[sweep2][scan2]['x'] = x
#            
#    del path,path2        
#
