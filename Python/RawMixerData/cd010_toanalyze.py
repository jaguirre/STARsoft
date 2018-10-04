# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 08:43:24 2018

@author: Alyssa
"""

from astropy import units as u

#cd010_dates = ['20180305',
#               '20180305',
#               '20180305',
#               '20180306',
#               '20180307',
#               '20180307',
#               '20180307',
#               '20180307',
#               '20180307',
#               '20180307',
#               '20180307',
#               '20180309',
#               '20180309',
#               '20180309',
#               '20180309',
#               '20180309',
#               '20180309',
#               '20180312',
#               '20180312',
#               '20180312',
#               '20180312',
#               '20180312',
#               '20180312']
#
#cd010_scans = ['noise13',
#               'noise15',
#               'noise17',
#               'noise0',
#               'noise0',
#               'noise01',
#               'noise02',
#               'noise03',
#               'noise04',
#               'noise05',
#               'noise06',
#               'noise0',
#               'noise03',
#               'noise04',
#               'noise05',
#               'noise06',
#               'noise07',
#               'noise0',
#               'noise02',
#               'noise04',
#               'noise06',
#               'noise10',
#               'noise12']
#
#cd010_Tstage = u.mK*[210,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,
#                215,]
#
#
#cd010_TBB = u.K*[6.4,
#             6.75,
#             7.1,
#             6.75,
#             7.1,
#             7.6,
#             8.1,
#             8.7,
#             9.3,
#             10.2,
#             11.1,
#             6.4,
#             6.4,
#             6.4,
#             6.4,
#             6.4,
#             6.4,
#             6.5,
#             6.75,
#             7.1,
#             7.6,
#             8.7,
#             9.3]
#%%
''' CD010 SHD mixer data '''
TBB_list = u.K*[6.65,
                7.1,
                7.45,
                7.8,
                8.45,
                9.0,
                9.85,
                10.6]

Tstage_list = u.mK*[.215,
                    .215,
                    .215,
                    .215,
                    .215,
                    .215,
                    .215,
                    .215]

dates_list = len(TBB_list)*['20180328']

scans_list = ['noise01',
              'noise02',
              'noise03',
              'noise04',
              'noise05',
              'noise06',
              'noise07',
              'noise09']

if (len(TBB_list)==len(Tstage_list)==len(dates_list)==len(scans_list))==False:
    print('Missed something! Input parameters are not the same length.')
#%%
from ReadMixerData import *

cool = 'CD010'
parentdatafolder = '/scr/starfire/labdata/'
parentoutfolder = '/scr/starfire/analysis/'

for ind in [7]: #for ind in np.arange(0,len(scans_list)):
    testdict = {}
    datafolder = parentdatafolder + dates_list[ind] + '/' + scans_list[ind] + '/'
    datescan = dates_list[ind] + '_' + scans_list[ind]
    outfolder = parentoutfolder + cool + '/' + datescan + '/'
    T_stage = Tstage_list[ind]
    T_BB = TBB_list[ind]
    #importmixerdata(testdict,T_stage,T_BB,cool,datafolder,outfolder,Pn=1,Fn=1,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5)
    
    importmixerfolder(testdict,T_stage,T_BB,cool,datafolder,outfolder,datescan,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5) 

    savedatadict(testdict,cool,datescan,outfolder)
