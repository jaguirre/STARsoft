# -*- coding: utf-8 -*-
"""
Created on Thu May 04 10:56:24 2017

@author: Alyssa
"""

import numpy as np
import matplotlib.pyplot as plt
import glob

# [f0	Qr	Qc	Qi	a_nl T_stage	T_BB]


filelist = glob.glob('201705_TempCharacterization/TempSweep/*.csv')

n_scans = 6
n_resonators = 12
n_parameters = 7

list = np.ndarray(shape=(n_scans,n_resonators,n_parameters))

fig = plt.figure(1) 
p1 = fig.add_subplot(311) # plot 1 -- f0 vs T
p2 = fig.add_subplot(312) # plot 2 -- Qr vs T
p3 = fig.add_subplot(313) # plot 3 -- Qi vs T

for scan in np.arange(0,n_scans):
    fit_results = np.loadtxt(filelist[scan],skiprows=2,delimiter=',')
    list[scan] = fit_results

colors = ['white','red','coral','orange','gold','y',
          'lightgreen','g','cyan','dodgerblue','darkorchid','hotpink']

# note: intentionally leaving out the first resonator here since the data are duplicated in the second
for res in np.arange(1,n_resonators):
    
    f0_list = list[:,res,0]
    temp_list = list[:,res,5]
    Qr_list = list[:,res,1]
    Qi_list = list[:,res,3]
    
    reslabel = 'Resonator ' + str(res)
  
    p1.plot(temp_list,f0_list,'o-',color=colors[res],label=reslabel)
    p1.set_xlabel('Stage Temperature (K)')
    p1.set_ylabel('Qr')
    
    p2.semilogy(temp_list,Qr_list,'o-',color=colors[res])
    p2.set_xlabel('Stage Temperature (K)')
    p2.set_ylabel('Qr')

    p3.semilogy(temp_list,Qi_list,'o-',color=colors[res])
    p3.set_xlabel('Stage Temperature (K)')
    p3.set_ylabel('Qi')
    
    
handles, labels = p1.get_legend_handles_labels()    
fig.legend(handles,labels)    
fig.show()
