#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:50:33 2017

@author: jaguirre
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c
import kid_responsivity_module as kids
reload(kids)
#%%

T = np.linspace(0.1,0.5,num=100)
#%%
resp1 = kids.kid_responsivity(T,0.1,vol=38.,alpha_k=0.2,nstar=100,tau_max=1e-3,nu_opt=850.)
resp2 = kids.kid_responsivity(T,0.1,vol=38.,alpha_k=0.2,nstar=1e4,tau_max=1e-5,nu_opt=850.)

#%%
plt.figure(1)
plt.clf()
plt.semilogy(T,resp1['nqp']+400.,'k')
plt.semilogy(T,resp2['nqp']+400.,'k--')
