#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 21:18:20 2018

@author: starfire
"""

from ReadMixerData import *
from AnalyzeResonatorData import *
import numpy as np

date,scan,T_stage,T_BB = np.loadtxt('CD010_toanalyze.csv',dtype=[str,str,float,float],skiprows=1,delimiter=',',unpack=True)



