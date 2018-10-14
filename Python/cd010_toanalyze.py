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
              'noise07']
              #'noise09']

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

cool = 'CD010'
parentdatafolder = '/scr/starfire/labdata/'
parentoutfolder = '/scr/starfire/analysis/'

#for ind in [7]: #
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
    
    
#%%
#%%
#plt.figure(100)
#colors = ['maroon','red','darkorange','green','limegreen','b','dodgerblue','m']
#
#res1 = 0.33163*u.GHz
#res2 = 0.34217*u.GHz
#res3 = 0.35553*u.GHz
#res4 = 0.36178*u.GHz
#resonators = [res1,res2,res3,res4]
#
#Qrange = 300
#df = [f0/Qrange for f0 in resonators]
#
#for res,resfreq in enumerate(resonators):
#    
#    res=0
#    resfreq=resonators[res]
#    fig,(p1,p2,p3) = plt.subplots(3,1,sharex=True,num=(str(res)+' LB'))
#    fig2,(px,psxx) = plt.subplots(2,1,sharex=True,num=(str(res)+' TBB'))
#    
#    for ind,scan in enumerate(scans_list[:-1]):
#    #ind=0
#    #scan=scans_list[ind]
#        testdict = hdf5_to_dict('/scr/starfire/analysis/CD010/20180328_'+scan+'/CD010_20180328_'+scan+'_datadict.hdf5')
#        
#        # Label each dataset by resonator
#        for pn in testdict.keys():
#            try:
#                for fn in testdict[pn].keys():
#    #                print(str(pn)+', '+str(fn))
#                    for indx,reso in enumerate(resonators):
#                        if abs(reso-testdict[pn][fn]['initial f0']*u.Hz) < df[indx]:
#                            testdict[pn][fn]['res'] = indx
#    
#            except AttributeError: exit
#        
#        ### want: P_LB vs a_nl,Qr,Qc,Qi,
#        # separate out by resonator
#    #    for res,resfreq in enumerate(resonators):
#        paths = dictwhere(testdict,'res',res)
#        streampaths = [[paths[k][0],paths[k][1],'stream'] for k in np.arange(len(paths))]
##        try:
#        p1.plot(dictget(testdict,paths,'LB_atten'),dictget(testdict,paths,'anl_fit'),'o',color=colors[ind],alpha=0.4)
#        p2.plot(dictget(testdict,paths,'LB_atten'),dictget(testdict,paths,'f0_fit'),'o',color=colors[ind],alpha=0.4)
#        p3.plot(dictget(testdict,paths,'LB_atten'),dictget(testdict,paths,'Qr_fit'),'s',color=colors[ind],alpha=0.4)
#        p3.plot(dictget(testdict,paths,'LB_atten'),dictget(testdict,paths,'Qc_fit'),'v',color=colors[ind],alpha=0.4)
#        
#        px.plot(dictget(testdict,paths,'T_BB'),dictget(testdict,paths,'f0_fit'),'ko',alpha=0.4)
#        psxx.plot(dictget(testdict,paths,'T_BB'),dictget(testdict,streampaths,'amp sub Sxx white'),'ko',alpha=0.4)
#    p1.set_xlim([25,45])
#    p1.set_ylim([-1,.25])
##    avals = list(gen_dict_extract('anl_fit',testdict))
##    pvals = list(gen_dict_extract('LB_atten',testdict))
##    plt.plot(pvals,avals,'.',color=colors[ind],alpha=0.5)
#    
##%%    
#    
#''' Resonance frequencies of resonators we studied for this device: '''
#res0 = 324.56*u.MHz
#res1 = 331.6*u.MHz
#res2 = 336.84*u.MHz
#res3 = 342.22*u.MHz
#res4 = 352.03*u.MHz
#res5 = 355.61*u.MHz
#res6 = 361.85*u.MHz
#res7 = 363.59*u.MHz
#resonators=[res0,res1,res2,res3,res4,res5,res6,res7]
#
#
#''' We'll use the range df=f0/Qrange to determine which resonator is which
#    [not the most elegant solution if there are lots of tightly-packed resonators and collisions] '''
#Qrange = 300
#df = [f0/Qrange for f0 in resonators]
#
#
#
#for cool in cdlist:
#    
#    ''' Compare each resonance frequency to the known list of resonators, 
#        use the value to index each scan by resonator number '''
#    for freq in list(gen_dict_extract('f0',cool)):
#        path = dictwhere(cool,'f0',freq)
#        for ind,res in enumerate(resonators):
#            if abs(res-freq) < df[ind]:
#                for sweep,scan in path: # hard-coding in that each freq will be indexed by cd[sweep][scan]
#                    cool[sweep][scan]['res'] = ind    
#    
#    ''' Find the maximum resonance frequency for each resonator within a cooldown,
#        set the maximum to be f00 for that resonator,cooldown '''
#    f00 = np.zeros(len(resonators))*u.MHz
#    for n,r in enumerate(resonators):
#        try: 
#            f0n = max(dictget(cool,dictwhere(cool,'res',n),'f0'))
#            f00[n] = f0n
#            
#        # if a resonator doesn't appear in this cooldown, set its f00 = -1
#        except ValueError: 
#            f00[n] = -1*u.MHz
#    cool['f00'] = f00
#    
#    ''' Use the f00 values for each resonator,cooldown to calculate x = (f-f00)/f
#        for each scan in the cooldown '''
#    for freq in list(gen_dict_extract('f0',cool)):
#        path2 = dictwhere(cool,'f0',freq)
#        reso = dictget(cool,path2,'res')
#        x = ((freq-f00[reso])/f00[reso]).to(u.dimensionless_unscaled)
#        for sweep2,scan2 in path2: # hard-coding in that each freq will be indexed by cd[sweep][scan]
#            cool[sweep2][scan2]['x'] = x
#            
#    del path,path2        
#
