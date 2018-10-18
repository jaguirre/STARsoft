#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 10:45:05 2018

@author: starfire
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import KID_model_functions as kids
import fitting_KID_model_functions as fitkids

res=0
for cool in ['CD010','CD011','CD012','CD013']:
    datafile = '/scr/starfire/analysis/'+cool+'/'+cool+'_reduced_Res'+str(res)+'.csv'
    
    TBBlist,f0list,xlist,Sxxlist = np.loadtxt(datafile,delimiter=',',unpack=True)
    
    
    