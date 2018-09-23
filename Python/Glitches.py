
# coding: utf-8

# In[98]:

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav
from astropy.stats import sigma_clip
import scipy.signal as sig


# In[2]:

good_datadir = '/scr/starfire/testdata/CD012/Analog/OpticalTests/mixer_data/20180607_noise03/Set0000Pn00Fn00/'
bad_datadir = '/scr/starfire/testdata/CD012/Analog/OpticalTests/mixer_data/20180607_noise03/Set0000Pn03Fn03/'
filename='calibrateddata.sav'


# In[3]:

dat_good = readsav(good_datadir+filename,python_dict=True)
dat_bad = readsav(bad_datadir+filename,python_dict=True)


# In[33]:

S21_good = dat_good['mixer_data_calibrated'][0]['S21_STREAM_CORR_ROT']
I_good = S21_good.real
Q_good = S21_good.imag
S21_bad = dat_bad['mixer_data_calibrated'][0]['S21_STREAM_CORR_ROT']
I_bad = S21_bad.real
Q_bad = S21_bad.imag


# In[34]:

t_good, dt_good = np.linspace(0,45,num=len(I_good),retstep=True)
t_bad, dt_bad = np.linspace(0,45,num=len(I_bad),retstep=True)


# In[62]:

I_bad_sc = sigma_clip(I_bad,sigma=5,iters=5)
Q_bad_sc = sigma_clip(Q_bad,sigma=5,iters=5)
good_data = ~np.logical_or(I_bad_sc.mask,Q_bad_sc.mask)


# In[63]:

good_data


# In[64]:

plt.figure(figsize=(20,8))
plt.plot(t_bad,I_bad)
plt.plot(t_bad[good_data],I_bad_sc.data[good_data],'o')
plt.plot(t_bad,Q_bad)
plt.plot(t_bad[good_data],Q_bad[good_data],'o')
plt.xlim([0.795,0.8])
plt.show()


# In[65]:

print(np.angle(S21_bad).min(),np.angle(S21_bad).max())
print(np.angle(S21_bad[~I_bad_sc.mask]).min(),np.angle(S21_bad[~I_bad_sc.mask]).max())


# In[122]:

p_bad,f_bad = plt.psd(S21_bad-np.mean(S21_bad),NFFT=2048,Fs=1/dt_bad)
p_good, f_good = plt.psd(S21_good-np.mean(S21_good),NFFT=2048,Fs=1/dt_good)
pI_good, fI_good = plt.psd(S21_good.real-np.mean(S21_good.real),NFFT=4096*4,Fs=1/dt_good)
pQ_good, fQ_good = plt.psd(S21_good.imag-np.mean(S21_good.imag),NFFT=4096*4,Fs=1/dt_good)
pI_bad, fI_bad = plt.psd(I_bad_sc-np.mean(I_bad_sc),NFFT=4096*4,Fs=1/dt_good)
pQ_bad, fQ_bad = plt.psd(Q_bad_sc-np.mean(Q_bad_sc),NFFT=4096*4,Fs=1/dt_good)
plt.psd(S21_good-np.mean(S21_good)-sig.convolve(S21_good,np.ones(11)/11.,mode='same'),NFFT=1024,Fs=1/dt_good,detrend='linear')
plt.show()


# In[124]:

plt.figure(figsize=(20,8))
plt.loglog(fI_good[fI_good>0],pI_good[fI_good>0],label='I good')
plt.loglog(fQ_good[fQ_good>0],pQ_good[fQ_good>0],label='Q good')
plt.loglog(fI_bad[fI_bad>0],pI_bad[fI_bad>0],label='I bad')
plt.loglog(fQ_bad[fQ_bad>0],pQ_bad[fQ_bad>0],label='Q bad')
plt.legend(loc='lower left')
plt.show()


# In[79]:

plt.semilogy(f_bad[f_bad>0],p_bad[f_bad>0])
plt.semilogy(f_good[f_good>0],p_good[f_good>0])
plt.show()


# In[87]:

plt.semilogy(np.fft.fftshift(np.abs(np.fft.fft(S21_good-np.mean(S21_good)))))
plt.ylim([1e-2,1e4])
plt.show()


# In[89]:

np.corrcoef(I_good,Q_good)


# In[102]:

plt.figure(figsize=(20,8))
plt.plot(t_good,I_good)
plt.plot(t_good,sig.convolve(S21_good,np.ones(11)/11.,mode='same').real)
plt.xlim([3,3.1])
plt.show()


# In[101]:

test = sig.convolve(S21_good,np.ones(11)/11.)


# In[106]:

get_ipython().magic('pinfo plt.psd')


# In[134]:

pI_bad_sc, fI_bad_sc = plt.psd(I_bad_sc-np.mean(I_bad_sc),NFFT=4096*4,Fs=1/dt_good)
pQ_bad_sc, fQ_bad_sc = plt.psd(Q_bad_sc-np.mean(Q_bad_sc),NFFT=4096*4,Fs=1/dt_good)
pI_bad, fI_bad = plt.psd(I_bad-np.mean(I_bad),NFFT=4096*4,Fs=1/dt_good)
pQ_bad, fQ_bad = plt.psd(Q_bad-np.mean(Q_bad),NFFT=4096*4,Fs=1/dt_good)
pI_bad_trim, fI_bad_trim = plt.psd((I_bad_sc.data-np.mean(I_bad_sc.data))[good_data],NFFT=4096*4,Fs=1/dt_good)


# In[142]:

plt.loglog(fI_bad_sc,pI_bad_sc.data)
plt.loglog(fI_bad,pI_bad)
plt.loglog(fI_bad_trim,pI_bad_trim)
plt.show()


# In[141]:

np.abs(pI_bad_trim-pI_bad_sc).max()


# /home/shd/IDL_Routines/my_routines/signal_processing/my_cross_spec.pro
# 
