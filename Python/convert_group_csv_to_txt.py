#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 23:36:03 2017

@author: jaguirre
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

#%%
# -------
# Things the user must set (should eventually come in from the command line or something)
csvfile = 'BBscan20170516.csv' # Name of your csv file.  Assumed to be in csvpath (defined below)
testdatapath = '/scr/starfire/testdata/'
CDlabel = 'CD004'
devlabel='BBscan20170516' # same devlabel as in process_one_group
filelabel='BBscan' # The prefix for the .txt files to be generatd.  date and time will be appended
# --------
# These are defined from above
csvpath=testdatapath+CDlabel+'/'
devpath = csvpath+devlabel+'/'

# Get the description of the groups from the csv file
descr = pd.read_csv(csvpath+csvfile)
try:
    os.system('mkdir '+devpath)
except:
    print('Cannot make devlabel directory')
    stop
os.system('mv '+csvpath+csvfile+' '+devpath)

#%%

def write_group_file(grp,CDlabel='CD004',devlabel='BBscan20170516',filelabel='',devpath=''):
    # Write the file that describes that group
    # grp is a pandas DataFrame
    
    # Generate the file name.  Name after stream.  Will fail if stream not present

    filename = filelabel+str(grp.loc[grp.type=='stream','date'].values[0])+'_'+ \
    str(grp.loc[grp.type=='stream','time'].values[0])+'.txt'
    #print filename

    f = open(devpath+filename,'w')
    f.write(CDlabel+'\n')
    f.write(devlabel+'\n')
    types = ['fine','gain','rough','med','stream']
    # The datelabels
    for t in types:
        #print t
        try:
            f.write(str(grp.loc[grp.type==t,'date'].values[0])+'\n')
        except:
            f.write('\n')
    # The timelabels
    for t in types:
        #print t
        try:
            f.write(str(grp.loc[grp.type==t,'time'].values[0])+'\n')
        except:
            f.write('\n')

    f.write('{0:0.004g}'.format(grp['T_bath'].median()/1000.)+'\n')
    f.write('{0:0.004g}'.format(grp['T_BB'].median())+'\n')

    f.close()

    return filename


#%%
for g in np.arange(1,descr['group'].max()+1):
    # Pick out a group 
    grp = descr.loc[descr.group==g] # pandas way of selecting a subset
    print(write_group_file(grp,CDlabel=CDlabel,devlabel=devlabel,filelabel=filelabel,devpath=devpath))
