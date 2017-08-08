#!/usr/bin/env python2                                                                                                                                                                                                             
# -*- coding: utf-8 -*-                                                                                                                                                                                                            
"""                                                                                                                                                                                                                                
Created on Thur June 29 10:45:50 2017                                                                                                                                                                                              
                                                                                                                                                                                                                                   
@author: Tasha                                                                                                                                                                                                                     
"""

import numpy as np
import astropy.units as u
from astropy.table import Table
import glob
import os
#from DataDictionary  import FinalDictionary                                                                                                                                                                                       
import pidly

idl = pidly.IDL()

fileconvert = 'COPYconvert_group_csv_to_txt.py'
datadictionary = 'DataDictionary.py'
plots = 'CreatePlots.py'
fileprocess = 'process_all_groups' #IDL Procedure                                                                                                                                                                                  

filelabel = 'InputPowerScan3'
#date = '20170725'                                                                                                                                                                                                                 
devlabel = filelabel#+date                                                                                                                                                                                                         
csvfile= devlabel+'.csv'

pathtodir = '/scr/starfire/testdata/CD006/'+ devlabel
#pathtofitresults = glob.glob(pathtodir+'/multitone/*/fit_results*.csv')                                                                                                                                                           
#pathtoCSV = 'cd /scr/starfire/testdata/'+ csvfile                                                                                                                                                                                 
#pathtoIDL = 'cd Tasha/STARsoft/IDL/test_pipelining/'                                                                                                                                                                              

print('Opening Python Script: '+ fileconvert + ' to begin converting csv file to text.')

"""                                                                                                                                                                                                                                
The goal of this section is to write a few lines of code that will replace the variables 'csvfile','devlabel', and 'filelabel' with the new names. Then you execute then csv-to-txt convert script.                                
"""

linenum = [16,19,20] # Line nmuber to be changed in csv-to-txt convert script.                                                                                                                                                     
linename = [csvfile,devlabel,filelabel]
linename2 = ['csvfile = ','devlabel = ','filelabel = ']

for num in np.arange(3):
    f = open(fileconvert,'r')
    lines = f.readlines()
    index = linenum[num]
    lines[index]= linename2[num] + "'" + linename[num] + "'\n" # Name of your csv file.  Assumed to be in csvpath                                                                                                                  
    f.close()
    f = open(fileconvert,'w')
    f.writelines(lines)
    f.close()
os.system('python '+fileconvert)

print('Finished converting csv file to txt file!')

"""                                                                                                                                                                                                                                
The goal of this section is to call the IDL procedure that processes all the groups of data taken.                                                                                                                                 
"""
idl('test')
print('Your Test worked so why not process_all_groups :(')

print('Opening IDL File: '+ fileprocess + ' to begin processing all scans.')
idl(fileprocess+ ',' + "'" + devlabel + "'")
print('Finished IDL execution!')

"""                                                                                                                                                                                                                                
The goal of this section is to write a few lines of code that will replace the variables 'path' and 'n_resonators'  when getting the fit*results.csv                                                                               
"""

f = open(datadictionary,'r')
lines = f.readlines()
lines[11] = 'path = '+ "'" + pathtodir + "'\n" # changes path address in DataDictionary.py                                                                                                                                         
f.close()
f = open(datadictionary,'w')
f.writelines(lines)
f.close()
os.system('python '+datadictionary)
print('Finished making Data Directory!')

"""                                                                                                                                                                                                                                
The goal of this section is to write a few lines of code that will replace the variables 'path' and 'n_resonators'  when getting the fit*results.csv                                                                               
"""
# I guess eventually we want to be able to change mutilple things for the plotting script...                                                                                                                                       

os.system('python '+ plots)
print('Finished making Plots!')