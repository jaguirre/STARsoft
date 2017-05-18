# -*- coding: utf-8 -*-
"""
Created on Thu May 04 10:56:24 2017

@author: Alyssa
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
#%%
# [f0	Qr	Qc	Qi	a_nl T_stage	T_BB]

path = '/scr/starfire/testdata/CD004/BBscan20170516/'
figfile = path+'BBscan20170516_fQ_vs_T.png'
filelist = glob.glob(path+'/multitone/*/fit_results*.csv')
#%%

n_scans = len(filelist)
n_resonators = 12
n_parameters = 7

data = np.ndarray(shape=(n_scans,n_resonators,n_parameters))

fig = plt.figure(1) 
plt.clf()
p1 = fig.add_subplot(311) # plot 1 -- f0 vs T
p2 = fig.add_subplot(312) # plot 2 -- Qr vs T
p3 = fig.add_subplot(313) # plot 3 -- Qi vs T

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
    Qi_list = data[:,res,3]
    
    indep_var = T_BB
    xlabel = 'BB Temperature (K)'
    
    reslabel = 'Resonator ' + str(res)
  
    p1.plot(indep_var,x,'o-',color=colors[res],label=reslabel)
    p1.set_xlabel(xlabel)
    p1.set_ylabel(r'$(f-f_0)/f_0$')
    
    p2.semilogy(indep_var,Qr_list,'o-',color=colors[res])
    p2.set_xlabel(xlabel)
    p2.set_ylabel('Qr')

    p3.semilogy(indep_var,Qi_list,'o-',color=colors[res])
    p3.set_xlabel(xlabel)
    p3.set_ylabel('Qi')
    
    
handles, labels = p1.get_legend_handles_labels()    
fig.legend(handles,labels)    
plt.savefig(figfile)
fig.show()
