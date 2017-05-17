#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 23:36:03 2017

@author: jaguirre
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%%

testfile='BBScan20170516.csv'

# Get the description of the groups
descr = pd.read_csv(testfile)

#%%

def write_group_file(grp,CDlabel='CD004',devlabel='BBscan20170516',filelabel=''):
    # Write the file that describes that group
    # grp is a pandas DataFrame
    
    # Generate the file name.  Name after stream.  Will fail if stream not present

    filename = filelabel+str(grp.loc[grp.type=='stream','date'].values[0])+'_'+ \
    str(grp.loc[grp.type=='stream','time'].values[0])+'.txt'
    #print filename

    f = open(filename,'w')
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
    print write_group_file(grp,filelabel='BBscan')