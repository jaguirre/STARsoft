# -*- coding: utf-8 -*-
"""
Created on Thu May 04 10:56:24 2017

@author: Alyssa
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy import units as u
from astropy import constants as c
#%%
# [f0	Qr	Qc	Qi	a_nl T_stage	T_BB]

#path = '/scr/starfire/testdata/CD004/BBscan20170516/'
path = '/Users/jaguirre/Documents/STARFIRE/BBscan20170516/'
figfile = 'fandQ_vs_T_BB.png'
filelist = glob.glob(path+'/multitone/*/fit_results*.csv')
#%%

n_scans = len(filelist)
n_resonators = 12
n_parameters = 7

data = np.ndarray(shape=(n_scans,n_resonators,n_parameters))

fig = plt.figure(1) 
plt.clf()
p1 = fig.add_subplot(411) # plot 1 -- f0 vs T
p2 = fig.add_subplot(412) # plot 2 -- Qr vs T
p3 = fig.add_subplot(413) # plot 3 -- Qi vs T
p4 = fig.add_subplot(414) # plot 4 -- Qc vs T; sanity check

for scan in np.arange(0,n_scans):
    fit_results = np.loadtxt(filelist[scan],skiprows=2,delimiter=',')
    data[scan] = fit_results

colors = ['white','red','coral','orange','gold','y',
          'lightgreen','g','cyan','dodgerblue','darkorchid','hotpink']

# note: intentionally leaving out the first resonator here since the data are duplicated in the second
for res in np.arange(1,n_resonators):
 
    f0_list = data[:,res,0]
    f0 = f0_list.max()
    
    x = (f0_list-f0)/f0
    T_stage = data[:,res,5]
    T_BB = data[:,res,6]
    Qr_list = data[:,res,1]
    Qc_list = data[:,res,2]
    Qi_list = data[:,res,3]
    
    indep_var = T_BB
    xlabel = 'BB Temperature (K)'
    
    reslabel = 'Resonator ' + str(res)
  
    p1.plot(indep_var,x,'o-',color=colors[res],label=reslabel)
    #p1.set_xlabel(xlabel)
    p1.set_ylabel(r'$(f-f_0)/f_0$')
    
    p2.semilogy(indep_var,Qr_list,'o-',color=colors[res])
    #p2.set_xlabel(xlabel)
    p2.set_ylabel(r'$Q_r$')

    p3.semilogy(indep_var,Qi_list,'o-',color=colors[res])
    #p3.set_xlabel(xlabel)
    p3.set_ylabel(r'$Q_i$')
    
    p4.semilogy(indep_var,Qc_list,'o-',color=colors[res])
    p4.set_xlabel(xlabel)
    p4.set_ylabel(r'$Q_c$')
      
handles, labels = p1.get_legend_handles_labels()    
fig.legend(handles,labels)    
plt.savefig(figfile)
fig.show()
#%%
psd_vals = np.loadtxt('fit_vals.csv',delimiter=',')
xall = data[:,1:,0].T
#%%
for res in np.arange(xall.shape[0]):
    f0 = xall[res,:].max()
    xall[res,:] = (xall[res,:]-f0)/f0 

#%%
plt.figure(2)
plt.clf()
Sxx=np.exp(psd_vals)
for res in np.arange(xall.shape[0]):
    resno = res+1 # fix the fact that Joe has already thrown out the first one
    reslabel = 'Resonator ' + str(resno)
    plt.plot(np.abs(xall[res,:]),Sxx[res,:],'o-',color=colors[resno],label=reslabel)

fit = np.polyfit(np.abs(xall).flatten(),Sxx.flatten(),1)
xp = np.linspace(0,0.0002)
bestfit = np.polyval(fit,xp)

#dS = 1e-16 * np.power(u.Hz,-1)
#dx = 2e-4
#m = dS/dx #1e-16/2e-4
#Sxx0 = 0.5e-16*np.power(u.Hz,-1)

dSdx = fit[0]*np.power(u.Hz,-1)
Sxx0 = fit[1]*np.power(u.Hz,-1)

nu = 857.*u.GHz
R = (dSdx/(4.*c.h*nu)).to(np.power(u.W,-1))
NEP = (np.sqrt(Sxx0)/R).to(u.W*np.power(u.Hz,-0.5))

plt.plot(xp,bestfit,'k',linewidth=2)
#plt.plot(xp,dSdx*xp+Sxx0,'k',linewidth=2)
plt.xlabel(r'|x| = $|\delta f/f|$')
plt.ylabel(r'$S_{xx}$')
plt.legend()
plt.savefig('Sxx_vs_x.png')
