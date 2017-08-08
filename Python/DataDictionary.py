# -*- coding: utf-8 -*-                                                                                                                                                                                                            
"""                                                                                                                                                                                                                                
Created on Wed July 05 1430:24 2017                                                                                                                                                                                                
@author: Tasha                                                                                                                                                                                                                     
"""

import numpy as np
from astropy.table import Table
import glob

# [f0   Qr      Qc      Qi      a_nl   T_stage  T_BB    atten]                                                                                                                                                                     
path = '/scr/starfire/testdata/CD006/InputPowerScan3'
filelist = glob.glob(path+'/multitone/*/fit_results*.csv') # We write to this line in another script                                                                                                                               
filelabel = 'InputPowerScan3' # We write to this line in another script                                                                                                                                                            

n_scans = len(filelist)
n_resonators = 11
n_parameters = 7 # This was 8 but atten is not included                                                                                                                                                                            
data = np.ndarray(shape=(n_scans,n_resonators,n_parameters))
indextoexclude = 'No' # Do you want to exlcude any indicies 'Yes' or 'No'                                                                                                                                                          

for scan in np.arange(0,n_scans):
    fit_results = np.loadtxt(filelist[scan],skiprows=2,delimiter=',')
    data[scan] = fit_results

array1_is_out_of_order = data[:,7] #Power sampling was not done in order so I'm fixing that here                                                                                                                                   
srt = np.argsort(array1_is_out_of_order)
array1_is_in_order = array1_is_out_of_order[srt]

"""                                                                                                                                                                                                                                
NOTE: intentionally leaving out the first resonator here since the data are duplicated in the second. I also wrote the next 3 lines to exclude any problem resonators.                                                            \
                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                   
If you know that there are resonators that you want to leave out of the plot then this if statement will be useful.                                                                                                                
"""
#if indextoexclude == 'Yes':                                                                                                                                                                                                       
#                                                                                                                                                                                                                                  
#    exclude_index = [8]                                                                                                                                                                                                           
#    a = np.arange(1,n_resonators)                                                                                                                                                                                                 
#    new_numRes = a[np.arange(len(a))!=exclude_index]                                                                                                                                                                              
#else :                                                                                                                                                                                                                            
#    a = np.arange(1,n_resonators)                                                                                                                                                                                                 
#    new_numRes = a                                                                                                                                                                                                                

FinalDictionary = {}

for res in np.arange(1,n_resonators): #new_numRes:                                                                                                                                                                                 

    f0_list = data[:,res,0].tolist()
    f0 = data[:,res,0].max()
    x = ((f0_list-f0)/f0).tolist()
    T_stage = data[:,res,5].tolist()
    T_BB = data[:,res,6].tolist()
    Qr_list = data[:,res,1].tolist()
    Qi_list = data[:,res,3].tolist()
    Qc_list = data[:,res,2].tolist()
    a_nl = data[:,res,4].tolist()
#    atten = data[:,res,7].tolist()                                                                                                                                                                                                

   # tables = Table([f0_list,x,T_stage,T_BB,Qr_list,Qi_list,Qc_list,a_nl,atten],names=('Freq List','Fractional Freq Shift','Stage Temp','BB Temp','Qr List','Qi List','Qc List','a_nl','Attenuation'),meta={'name':'first table'})
    tables = Table([f0_list,x,T_stage,T_BB,Qr_list,Qi_list,Qc_list,a_nl],names=('Freq List','Fractional Freq Shift','Stage Temp','BB Temp','Qr List','Qi List','Qc List','a_nl'),meta={'name':'first table'})
    FinalDictionary['Res '+ str(res)] = tables
#print(FinalDictionary)        