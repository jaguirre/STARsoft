# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 08:43:24 2018

@author: Alyssa
"""

from astropy import units as u

#%%
''' CD011 SHD mixer data '''

Tstage_list = u.mK*[.215,
                    .215,
                    .225,
                    .25,
                    .275,
                    .3,
                    .325,
                    .35]

TBB_list = len(Tstage_list)*[u.K*0]

dates_list = len(Tstage_list)*['20180509']

scans_list = ['noise04',
              'noise05',
              'noise06',
              'noise07',
              'noise08',
              'noise09',
              'noise10']
#              'noise11']

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
#cool = 'CD011'
#parentdatafolder = '/scr/starfire/labdata/'
#parentoutfolder = '/scr/starfire/analysis/'
#
#coolfolder = parentoutfolder+cool+'/'
#if os.path.isdir(coolfolder) == False: os.mkdir(coolfolder)
#
#for ind in np.arange(0,len(scans_list)):
#    
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
trans=1
eta_opt = 0.17*u.dimensionless_unscaled
N0=1.72e10*np.power(u.micron,-3)*np.power(u.eV,-1)



n = len(scans_list)

for indx,res in enumerate(resonators):
    res=resonators[indx]
    Tstagelist = []
    f0list = []
    Qrinvlist = []
    Sxxlist = []
    alist = []
    LBlist = []

    for ind in np.arange(0,n):
        scan = scans_list[ind]
        testdict = hdf5_to_dict('/scr/starfire/analysis/CD011/20180509_'+scan+'/CD011_20180509_'+scan+'_datadict.hdf5')

        allpaths = dictwhere(testdict,'res',indx)    
        a = dictget(testdict,allpaths,'LB_atten')
        
        paths = [allpaths[pa] for pa in np.where(a==np.max(a))[0]]
        streampaths = [[paths[k][0],paths[k][1],'stream'] for k in np.arange(len(paths))]
        
        try:
            T_stage = dictget(testdict,paths,'T_stage')
            for T in T_stage: Tstagelist.append(T)
            
            f0 = dictget(testdict,paths,'f0_fit')
            for f in f0: f0list.append(f)
            
            Qr = dictget(testdict,paths,'Qr_fit')
            for Q in Qr: Qrinvlist.append(1./Q)
            
            Sxx = dictget(testdict,streampaths,'amp sub Sxx white')
            for S in Sxx: Sxxlist.append(S)
            
            a1 = dictget(testdict,paths,'anl_fit')
            for anl in a1: alist.append(anl)
            
            plb = dictget(testdict,paths,'LB_atten')
            for polb in plb: LBlist.append(polb)


        except KeyError: exit
        
    xlist = ((f0list-np.max(f0list))/np.max(f0list))
    Qrinvlist = np.array(Qrinvlist)
    Tstagelist = u.K*np.array(Tstagelist)
    xtest = kids.xMB(alpha,u.Hz*f0list[0],Tstagelist,Tc,TBB,V,n_star,tau_max,eta_pb,eta_opt,trans=0,N0=N0)    
    Qinvtest = kids.QinvMB(alpha,u.Hz*f0list[0],Tstagelist,Tc,6*u.K,V,n_star,tau_max,eta_pb,eta_opt=1,trans=0,N0=N0)
    data = [Tstagelist,f0list[0]*u.Hz]
    xsig = 0.05*xlist[0]*np.ones_like(xlist)
    if xlist[0]==0: xsig = 0.05*xlist[1]*np.ones_like(xlist)
    Qsig = 0.05*Qrinvlist
    
    xfitopt,xfitcov = curve_fit(fitkids.x_dark_fit,data,xlist,p0=(.73,1.39,(xtest[0]/(f0list[0]))),sigma=xsig)#,bounds=([0,.5,-1],[1,2.5,0]))
    Qfitopt,Qfitcov = curve_fit(fitkids.Qinv_dark_fit,data,Qrinvlist,p0=(.73,1.39,Qinvtest[0]),sigma=Qsig,bounds=([0,.5,0],[1,2.5,max(Qrinvlist)]))
    
    cdata = np.concatenate((xlist,Qrinvlist))
    csigma = np.concatenate((xsig,Qsig))
    p0=(.73,1.39,xfitopt[-1],Qfitopt[-1])#(xtest[0]/f0list[0]),Qinvtest[0])
    cfitopt,cfitcov = curve_fit(fitkids.x_Qinv_dark_simulfit,data,cdata,sigma=csigma,p0=(.51,1.32,4.1e-14,8.9e-6),bounds=([0,.5,-1,0],[1,2.5,1,max(Qrinvlist)]))
    print(str(res)+': ')
    print(xfitopt)
    print(Qfitopt)
    print(cfitopt)

    fig,(p1,p2,p3,p4) = plt.subplots(4,1,sharex=True,num=('Res '+str(indx)))
    fig.suptitle('Res '+str(indx))
    p1.plot(Tstagelist,xlist,'ko')
    p1.plot(Tstagelist,xtest,'g--')
    p1.plot(Tstagelist,fitkids.x_dark_fit(data,*xfitopt),'cx')
    p1.plot(Tstagelist,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[0:len(Tstagelist)],'r-')
    p1.set_ylabel(r'$\delta$ f/f')
    p2.plot(Tstagelist,Qrinvlist,'ko')
    p2.plot(Tstagelist,Qinvtest,'g--')
    p2.plot(Tstagelist,fitkids.Qinv_dark_fit(data,*Qfitopt),'cx')
    p2.plot(Tstagelist,fitkids.x_Qinv_dark_simulfit(data,*cfitopt)[len(Tstagelist):],'r-')
    p2.set_ylabel(r'1/$Q_r$')
    p2.set_title(r'Simul fit: $\alpha$ = {:.2f}'.format(cfitopt[0])+r', $T_c$ = {:.2f}'.format(cfitopt[1]))
    p3.plot(Tstagelist,Sxxlist,'ko')
    p3.set_ylabel(r'$S_{xx}$ (Hz$^{-1}$)')
    p4.set_xlabel(r'$T_{stage}$ (K)')
    p4.plot(Tstagelist,alist,'ko')
    p4.set_ylabel(r'$a_{nl}$ fit')
    p4b = p4.twinx()
    p4b.plot(Tstagelist,LBlist,'gs')
    p4b.set_ylabel(r'LB atten (dB)')
    
    