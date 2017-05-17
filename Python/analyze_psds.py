#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 16:09:16 2017

@author: jaguirre
"""

import psd_utils as psd


temps = [252, 275, 300, 325, 350, 375]
temp_dict = {temp:'{0}_tsweep_psddata.sav'.format(temp) for temp in temps}
psd.temp_data(temp_dict)