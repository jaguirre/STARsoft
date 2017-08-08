# -*- coding: utf-8 -*-                                                                                                                                                                                                            
"""                                                                                                                                                                                                                                
Created on Tue July 05 1559:15 2017                                                                                                                                                                                                
@author: Tasha                                                                                                                                                                                                                     
"""

import numpy as np
import matplotlib.pyplot as plt
from DataDictionary import FinalDictionary
import glob

# note: intentionally leaving out the first resonator here since the data are duplicated in the second. I also wrote the next 3 lines to exclude any problem resonators.                                                           

n_resonators = 11
indextoexclude = 'No' # Are there any resonators you want to exclude 'Yes' or 'No'.                                                                                                                                                
pathtodir = '/scr/starfire/testdata/CD006/InputPowerScan3/'

#if indextoexclude == 'Yes':                                                                                                                                                                                                       
#                                                                                                                                                                                                                                  
#    exclude_index = [8]                                                                                                                                                                                                           
#    a = np.arange(1,n_resonators)                                                                                                                                                                                                 
#    new_numRes = a[np.arange(len(a))!=exclude_index]                                                                                                                                                                              
#else :                                                                                                                                                                                                                            
#    a = np.arange(1,n_resonators)                                                                                                                                                                                                 
#    new_numRes = a                                                                                                                                                                                                                

d = FinalDictionary

fig = plt.figure(1)
fig.set_size_inches(12.5, 8.5,forward=True)#Resize figures in matplotlib                                                                                                                                                           
plt.clf()
p1 = fig.add_subplot(331) # plot 1 -- f0 vs indep_var                                                                                                                                                                              
p2 = fig.add_subplot(332) # plot 2 -- Qr vs indep_var                                                                                                                                                                              
p3 = fig.add_subplot(333) # plot 3 -- Qi vs indep_var                                                                                                                                                                              
p4 = fig.add_subplot(334) # plot 4 -- Qc vs indep_var                                                                                                                                                                              
p5 = fig.add_subplot(335) # plot 5 -- a_nl vs indep_var                                                                                                                                                                            

colors = ['white','red','coral','orange','gold','y',
          'lightgreen','g','cyan','dodgerblue','darkorchid','hotpink']

for res in np.arange(1,n_resonators):#new_numRes:                                                                                                                                                                                  
#    for s in np.arange(numScans):                                                                                                                                                                                                 

    indep_var = d['Res '+str(res)]['BB Temp']
    xlabel = 'BB Temperature (K)'# 'Attenuation (dB)'                                                                                                                                                                              
    reslabel = 'Resonator ' + str(res)

    p1.plot(indep_var,d['Res '+str(res)]['Fractional Freq Shift'],'-o',color=colors[res],label=reslabel)
    p1.set_xlabel(xlabel)
    p1.set_ylabel(r'$(f-f_0)/f_0$')

    p2.semilogy(indep_var,d['Res '+str(res)]['Qr List'],'-o',color=colors[res])
    p2.set_xlabel(xlabel)
    p2.set_ylabel('Qr')

    p3.semilogy(indep_var,d['Res '+str(res)]['Qc List'],'-o',color=colors[res])
    p3.set_xlabel(xlabel)
    p3.set_ylabel('Qc')

    p4.semilogy(indep_var,d['Res '+str(res)]['Qi List'],'-o',color=colors[res])
    p4.set_xlabel(xlabel)
    p4.set_ylabel('Qi')

#    p5.plot(indep_var,a_nl,'o',color=colors[res])                                                                                                                                                                                 
    p5.semilogy(indep_var,d['Res '+str(res)]['a_nl'],'-o',color=colors[res])
    p5.set_xlabel(xlabel)
    p5.set_ylabel('a_nl')

handles, labels = p1.get_legend_handles_labels()
fig.legend(handles,labels,bbox_to_anchor=(.95,.6),borderaxespad=0.) #,loc='best')                                                                                                                                                  
#fig.legend(bbox_to_anchor=(1.05, 1),handles,labels,loc=2, borderaxespad=0.)                                                                                                                                                       
#fig.subplots_adjust(hspace=.4,left=.2,bottom=None, right=None, top=None, wspace=None,)                                                                                                                                            
fig.tight_layout()
plt.savefig(pathtodir+'Plots.png')
#fig.show()             