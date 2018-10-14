# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 08:43:24 2018

@author: Alyssa
"""

from astropy import units as u

#%%
''' CD011 SHD mixer data '''


TBB_list = u.K*[5.61,
                5.75,
                5.95,
                6.65,
                7.1,
                7.45,
                7.8]

Tstage_list = len(TBB_list)*[215*u.mK]


dates_list = len(TBB_list)*['20180607']

scans_list = ['noise01',
              'noise02',
              'noise03',
              'noise04',
              'noise05',
              'noise06',
              'noise07']

res1 = 0.33163*u.GHz
res2 = 0.34217*u.GHz
res3 = 0.35553*u.GHz
res4 = 0.36178*u.GHz
resonators = [res1,res2,res3,res4]

Qrange = 300
df = [f0/Qrange for f0 in resonators]


if (len(TBB_list)==len(Tstage_list)==len(dates_list)==len(scans_list))==False:
    print('Missed something! Input parameters are not the same length.')
#%%
from ReadMixerData import *

cool = 'CD012'
parentdatafolder = '/scr/starfire/labdata/'
parentoutfolder = '/scr/starfire/analysis/'

coolfolder = parentoutfolder+cool+'/'
if os.path.isdir(coolfolder) == False: os.mkdir(coolfolder)

for ind in np.arange(0,len(scans_list)):
    
    testdict = {}
    datafolder = parentdatafolder + dates_list[ind] + '/' + scans_list[ind] + '/'
    datescan = dates_list[ind] + '_' + scans_list[ind]
    outfolder = parentoutfolder + cool + '/' + datescan + '/'
    T_stage = Tstage_list[ind]
    T_BB = TBB_list[ind]
    #importmixerdata(testdict,T_stage,T_BB,cool,datafolder,outfolder,Pn=1,Fn=1,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5)
    
    importmixerfolder(testdict,T_stage,T_BB,cool,datafolder,outfolder,datescan,docal=True,doPSD=True,doplots=True,Qignore=10**3,poly_order=5) 

    # label the resonators
    for pn in testdict.keys():
        try:
            for fn in testdict[pn].keys():
                for indx,reso in enumerate(resonators):
                    if abs(reso-testdict[pn][fn]['initial f0']*u.Hz) < df[indx]:
                        testdict[pn][fn]['res'] = indx

        except AttributeError: exit


    savedatadict(testdict,cool,datescan,outfolder)
