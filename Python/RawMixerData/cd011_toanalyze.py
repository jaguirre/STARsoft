# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 08:43:24 2018

@author: Alyssa
"""

from astropy import units as u

#%%
''' CD011 SHD mixer data '''

Tstage_list = u.mK*[.215,
                    .215,
                    .225,
                    .25,
                    .275,
                    .3,
                    .325,
                    .35]

TBB_list = len(Tstage_list)*[u.K*0]

dates_list = len(Tstage)*['20180509']

scans_list = ['noise04',
              'noise05',
              'noise06',
              'noise07',
              'noise08',
              'noise09',
              'noise10',
              'noise11']

if (len(TBB_list)==len(Tstage_list)==len(dates_list)==len(scans_list))==False:
    print('Missed something! Input parameters are not the same length.')
#%%
from ReadMixerData import *

cool = 'CD010'
parentdatafolder = '/scr/starfire/labdata/'
parentoutfolder = '/scr/starfire/analysis/'

for ind in np.arange(0,len(scans_list)):
    
    testdict = {}
    datafolder = parentdatafolder + dates_list[ind] + '/' + scans_list[ind] + '/'
    datescan = dates_list[ind] + '_' + scans_list[ind]
    outfolder = parentoutfolder + cool + '/' + datescan + '/'
    T_stage = Tstage_list[ind]
    T_BB = TBB_list[ind]
    #importmixerdata(testdict,T_stage,T_BB,cool,datafolder,outfolder,Pn=1,Fn=1,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5)
    
    importmixerfolder(testdict,T_stage,T_BB,cool,datafolder,outfolder,datescan,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5) 

    savedatadict(testdict,cool,datescan,outfolder)
