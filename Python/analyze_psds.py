#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 16:09:16 2017

@author: jaguirre
"""

import psd_utils as psd
from glob import glob
from scipy.io import readsav
import numpy as np

path = '/Users/jaguirre/Documents/STARFIRE/BBscan20170516/'
psdfilelist = glob(path+'/multitone/*/psddata*.sav')

#%%

# From the summary csv file
temps = data[:,0,6] #[252, 275, 300, 325, 350, 375]
#temp_dict = {temp:'{0}_tsweep_psddata.sav'.format(temp) for temp in temps}
temp_dict = {temp:psdfilelist[i] for i,temp in enumerate(temps)}
psd.temp_data(temp_dict)

#psddata = psd.read_data(temp_dict)
#data_sxx = [(key, psd.get_sxx(temp_data)) for key, temp_data in psddata.items()]

#%%
fit_vals = np.loadtxt('fit_vals.csv',delimiter=',')