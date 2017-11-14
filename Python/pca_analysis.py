#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 14:29:49 2017

@author: tashaleebillings
"""

#+
# NAME:
#       PCA_ANALYSIS
#
# PURPOSE:
#       Given a collection of datasets, execute a principle component analysis. This closely follows the PCA example from Coyote:
#         http://www.idlcoyote.com/code_tips/pca.html
#         
#       The input is a 2D array in which each column is a different dataset. The first step is to subtract the mean from each 
#       column, and optionally scale by the stddev if the /normalize_stdev keyword is set.
#       
#       The second step is to compute the covariance matrix, and the associated eigenvectors and eigenvalues. A check on the 
#       eigenvector/eigenvalue computation is executed. The data is then rotated into the basis set defined by the eigenvectors;
#       the important point here is that the covariance matrix in this new basis is diagonal.
#       
#       We then zero the principle components with the largest eigenvalues (or magnitude of eigenvalues if the /absolute keyword
#       is set), with the number of removed components set by the /nrem keyword. We then de-rotated this filtered data back into
#       the original frame.
#       
#       We check the eigenvectors, eigenvalues, and rotated data against the same computation done with the PCOMP.pro routine.
#       We do the same for PCA.pro, with the extra comparison of the covariance matrix. Note that we do the individual steps here,
#       rather than relying on the PCOMP.pro or PCA.pro routine, such that the rotation and de-rotation of the data is more transparent.
#
# CATEGORY:
#
# CALLING SEQUENCE:
#
# INPUTS:
#       data -- the m x n array to be analyzed; each of m columns is a dataset with n samples
#
# OPTIONAL INPUTS:
#
# KEYWORD PARAMETERS:
#       normalize_stdev -- if specified, scale each column of the input data array so as to have unit variance
#       absolute -- if specified, the eigenvalues (and corresponding eigenvectors) are sorted in descending order by their absolute value, rather than by their signed value {default is 0}
#       nrem -- the number of principle components to remove (those with the largest eigenvalues are removed) {default is 0}
#       tol -- the tolerance for comparing numbers to 0 {default is 1.e-12}
#       docompare -- if set, compare key results here with PCOMP.pro and PCA.pro
#
# OUTPUTS:
#       struct -- a structure with various computed quantities
#       ['initialdata'] -- the data analyzed here; equal to the input data, but with each column mean-subtracted, and optionally scaled to have unit variance 
#       ['covmatrix'] -- the covariance matrix of initialdata
#       ['evals'] -- eigenvalues of the covariance matrix
#       ['evecs'] -- eigenvectors of the covariance matrix (each ROW is a different eigenvector), normalized to unit magnitude
#       ['rotateddata'] -- initialdata rotated into the frame defined by the principle components
#       ['filtereddata'] -- initialdata with the specified number of principle components removed
#       ['filtereddata_restored'] -- filtereddata with the mean and variance of the initial data restored
#
# RESTRICTIONS:
#       x) I don't know when we might have a negative eigenvalue, so it's not clear how/when the /absolute keyword should be used
#
# EXAMPLE:
#
# MODIFICATION HISTORY:  
#
#
#-
#***********************************************************************
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from attrdict import AttrDict
from scipy import signal


def pca_analysis(data, normalize_stdev=1,absolute=0,nrem=1,tol=1.e-12):

# Define the new variable initialdata, equal to the initial data
# - we introduce a new variable so we can modify it in this function without changing the variable in the main program
# Compute the size of initialdata, and explicitly define m and n 
#IMPORTANT: You need to look at your data to see how it is behaving because you might need cut off some of the bad data at the end.
    #initialdata = data
    initialdata = data[0:28600,:] #This line removes the junk we see at the end of the data
    sz = initialdata.shape
    m = sz[1] # Number of Columns
    n = sz[0] # Number of Rows

# Check that nrem is less than (m-1)
    if nrem > (m-1):
        print('nrem =', nrem)
        print('m =', m)
        print('Cannot remove more principle components than number of independent data sets.')

# Subtract the mean from each cloumn of initialdata
# Optionally scale initialdata wuch that each column has unit variance
    mean_vec = np.mean(initialdata, axis=0) #Takes the mean of elements in same column. So m stays the same. Set axis=1 to get mean of each column vector.
    initialdata = initialdata - np.resize(mean_vec,(n,m)) # len(mean_vec) should be the number of columns, m.

    if normalize_stdev is True:
        norm_vec = np.std(initialdata,axis=0,ddof=1) #Takes the std of all elements in the same column
        initialdata = initialdata/np.resize(norm_vec,(n,m))

# Compute the covariance matrix
# This is an m x m array, with element (i,j) equal to the covariance of columns i and j in inputdata
# Default rowvar = True, if so initialdata.transpose() must be done if you set rowvar = False the column is the variable

    covmatrix = np.cov(initialdata, rowvar = False) #This produces an n x n matrix. Each column vector must be the same length.
    corrmatrix = np.corrcoef(initialdata,rowvar = False)

# Compute the eigenvalues and eigenvectors of the covariance matrix
# The result is a m-element array of eigenvalues
# The variable set by the /eigenvectors keyword is an m x m array, with each ROW a different eigenvector
#  - recall these eigenvectors are normalized to unit length
# The variable set by the /residual keyword is an m x m array, equal to (A*x-lambda*x), where:
#  - A is the covariance matrix
#  - x is the eigenvector
#  - lambda is the eigenvalue
#  - this array should be zero, and can be used as a check on the calculation
# The /absolute keyword specifies that the eigenvalues (and corresponding eigenvectors) are sorted by decreasing absolute value of the eigenvalue
    evals,evecs  = np.linalg.eig(covmatrix)

# Check the computation of the eigenvalues through the residual keyword


# Use the eigenvectors to transform the data into the principle components
#  - Recall this is a rotation, and the principle components are such that the rotated datasets are independent
#  - in other words, the covariance matrix of the rotated data is diagonal
# Note that I have transposed the following two equations relative to the Coyote article, but it is consistent
    featurevector = evecs.transpose()
    rotateddata = (np.inner(featurevector,initialdata)).transpose() # This is the matrix with the eigenvectors in the columns transposed so that the eigenvectors are now in the rows.

# Define a filter that removes the principle components with the largest eigenvalues
    filtermatrix = np.ones((n,m))
    #if nrem > 0:
    #    filter[:,0:nrem] = 0

    largest = np.where(evals == evals.max())[0]
    filtermatrix[:,largest] = 0.

# Multiple the rotated data by the filter, and derotate back into the original frame
    fv_inv = np.linalg.inv(featurevector)
    filtereddata = (np.inner(fv_inv,np.multiply(rotateddata,filtermatrix))).transpose()

# Scale and shift the filtered data to restore the mean and variance of the initial dataset
    filtereddata_restored = filtereddata
    if normalize_stdev is True:
        filtereddata_restored = filtereddata_restored*np.resize(norm_vec,(n,m))

    filtereddata_restored = filtereddata_restored# + np.resize(mean_vec,(n,m))
    
    struct = {'initialdata':initialdata, 'covmatrix':covmatrix, 'corrmatrix':corrmatrix, 'evals':evals, 'evecs':evecs,
              'rotateddata':rotateddata, 'filtereddata':filtereddata, 'filtereddata_restored':filtereddata_restored, 'filter':filtermatrix}
          
    return struct

def isblind(filename):
     savdata = readsav(filename,python_dict=True)
     allbins = savdata['multitone_data_raw'][0]['bins']
     blindbins = savdata['multitone_data_raw'][0]['blindbin']
 
     isblind = np.zeros_like(allbins)
 
     for m in allbins:
         for n in blindbins:
             if m == n:
                 index = np.where(allbins==m)
                 isblind[index] = 1
     return isblind

def findresonators(filename):

    data = readsav(filename,python_dict=True)
    freqs = data['multitone_data_raw'][0]['frlist']

    blinds = isblind(filename)

    resonatorlist = np.zeros(len(blinds)-np.sum(blinds))

    res = 0
    for ind in np.arange(0,len(blinds)):
        if blinds[ind] == 0:
            resonatorlist[res] = freqs[ind]
            res+=1

    return {'reslist': resonatorlist, 'isblind': blinds}

def reverse_numeric(x,y):
    return y-x

#%%
rawsav = ['/Users/tashaleebillings/Desktop/morereadingmaterial/InputPowerScanRawData.sav','/Users/tashaleebillings/Desktop/morereadingmaterial/HighPower1HighPassRawData.sav','/Users/tashaleebillings/Desktop/morereadingmaterial/HighPower2RawData.sav','/Users/tashaleebillings/Desktop/morereadingmaterial/LowPowerRawData.sav']
rawsavdir = {'InputPowerScanRawData':'/Users/tashaleebillings/Desktop/morereadingmaterial/InputPowerScanRawData.sav','HighPower1HighPassRawData' :'/Users/tashaleebillings/Desktop/morereadingmaterial/HighPower1HighPassRawData.sav','HighPower2RawData':'/Users/tashaleebillings/Desktop/morereadingmaterial/HighPower2RawData.sav','LowPowerRawData':'/Users/tashaleebillings/Desktop/morereadingmaterial/LowPowerRawData.sav'}
dataname = ['InputPowerScanRawData','HighPower1HighPassRawData','HighPower2RawData','LowPowerRawData']

i=2 # i goes from 0 to len(rawsav). I'm sure you can use a for loop to iterate through the .sav files but I did it one by one.
rawdata = readsav(rawsav[i])
blind = isblind(rawsav[i])
resFreq = findresonators(rawsav[i])
    
#rawdata = readsav(rawsav[i]) #/scr/starfire/testdata/CD006/InputPowerScan/0_multitone/20170719_183039/rawdata.sav')
#blind = isblind('/Users/tashaleebillings/Desktop/morereadingmaterial/InputPowerScanRawData.sav') #/scr/starfire/testdata/CD006/InputPowerScan/0_multitone/20170719_183039/rawdata.sav')

#If you look at the data you will see that there are 'nan' at the end so I set them equal to zero but if you look at line 82 for the variable 'initialdata' I cut those off because they give a large spike. I guess you can write code that goes in and identifies and removes the 'nan'.
I = (rawdata.multitone_data_raw.streamdata[0].stream_data_concat[0].s21i[0])
Q = (rawdata.multitone_data_raw.streamdata[0].stream_data_concat[0].s21q[0])
hm = np.min(np.where(np.logical_not(np.isfinite(Q))))
#I = I[:hm,:].copy()
#Q = Q[:hm,:].copy()
I[np.logical_not(np.isfinite(Q))] = 0. #Beacause I already cute the data off in line 82 technically I don't need this?
Q[np.logical_not(np.isfinite(Q))] = 0. #Beacause I already cute the data off in line 82 technically I don't need this?

Ires = np.zeros((len(np.isfinite(I)),len(blind)-np.sum(blind)))
Qres = np.zeros((len(np.isfinite(Q)),len(blind)-np.sum(blind)))
Iblind = np.zeros((len(np.isfinite(I)),np.sum(blind)))
Qblind = np.zeros((len(np.isfinite(Q)),np.sum(blind)))
I_newElem = np.zeros(len(blind)-np.sum(blind))
row = 0
w = np.where(blind==0)[0]
wb = np.where(blind==1)[0]
 
for row in np.arange(0,len(Ires)):
    Ires[row] = I[row][w]
    Qres[row] = Q[row][w]
    Iblind[row] = I[row][wb]
    Qblind[row] = Q[row][wb]
    row += 1

#Ires = I[:,w].copy()
#Iblind = I[:,wb].copy()
#Qres = I[:,w].copy()
#Qblind = I[:,wb].copy()

IQres = np.concatenate((Ires,Qres),axis=1)
IQblind = np.concatenate((Iblind,Qblind),axis=1)

pcaI = pca_analysis(data=Ires, normalize_stdev=1,absolute=0,nrem=1,tol=1.e-12)
pcaQ = pca_analysis(data=Qres, normalize_stdev=1,absolute=0,nrem=1,tol=1.e-12)
pcaIQ = pca_analysis(data=IQres, normalize_stdev=1,absolute=0,nrem=1,tol=1.e-12)

pcaIblind = pca_analysis(data=Iblind, normalize_stdev=1,absolute=0,nrem=1,tol=1.e-12)
pcaQblind = pca_analysis(data=Qblind, normalize_stdev=1,absolute=0,nrem=1,tol=1.e-12)
pcaIQblind = pca_analysis(data=IQblind, normalize_stdev=1,absolute=0,nrem=1,tol=1.e-12)

#Covariance is just an unstandardized version of correlation.  To compute any correlation, we divide the covariance by the standard deviation of both variables to remove units of measurement.
#So a covariance is just a correlation measured in the units of the original variables.

#print(pcaI['evecs'].shape)
#print(pcaI['covmatrix'].shape)
#print(pcaI['corrmatrix']) #constrained to being between -1 and 1 to assess the strength of a realtionship between two variables
#print(pcaI['covmatrix']) 

#%%
# MAKE PLOTS

plt.figure(i)
plt.clf()
plt.title('Res Evals:'+dataname[i])
plt.plot(sorted(pcaIQ['evals'],reverse=True), '.-')#,cmp=reverse_numeric))
plt.tight_layout()

fig, axes = plt.subplots(nrows=2 , ncols=2, figsize=(14,14))

c1 = axes[0,0].imshow(pcaIQblind['corrmatrix'])
fig.colorbar(c1,ax=axes[0][0])
axes[0,0].set_title('IQ Blind Tones Correlation')

c2 = axes[1,0].imshow(pcaIQblind['covmatrix'])
fig.colorbar(c2,ax=axes[1][0])
axes[1,0].set_title('IQ Blind Tones Covariance')

c3 = axes[0,1].imshow(pcaIQ['corrmatrix'])
fig.colorbar(c3,ax=axes[0][1])
axes[0,1].set_title('IQ Tones Correlation')

c4 = axes[1,1].imshow(pcaIQ['covmatrix'])
fig.colorbar(c4,ax=axes[1][1])
axes[1,1].set_title('IQ Tones Covariance')

fig.suptitle(dataname[i])
#handles, labels = axes.get_legend_handles_labels()
#fig.legend(handles,labels,bbox_to_anchor=(.95,.6),borderaxespad=0.)
#fig.tight_layout()

fig,axes = plt.subplots(nrows=3, ncols=2, figsize=(14,14))

p1 = axes[0,0].plot(pcaIQ['initialdata'][:,0],pcaIQ['initialdata'][:,11],'.')
axes[0,0].set_title('Initial Data')
axes[0,0].axis('equal')
#axes[0,0].set_xlabel('I')
#axes[0,0].set_ylabel('Q')

p2 = axes[1,0].plot(pcaIQ['initialdata'][:,1],pcaIQ['initialdata'][:,12],'.')
axes[1,0].axis('equal')
#axes[1,0].set_title('Res 2 Initial Data')
#axes[1,0].set_xlabel('I')
#axes[1,0].set_ylabel('Q')

p3 = axes[2,0].plot(pcaIQ['initialdata'][:,2],pcaIQ['initialdata'][:,13],'.')
axes[2,0].axis('equal')
#axes[2,0].set_title('Res 3 Initial Data')
#axes[2,0].set_xlabel('I')
#axes[2,0].set_ylabel('Q')

p4 = axes[0,1].plot(pcaIQ['filtereddata_restored'][:,0],pcaIQ['filtereddata_restored'][:,11],'.')
axes[0,1].set_title('Filtered Data Restored')
axes[0,1].axis('equal')
#axes[0,1].set_xlabel('I')
#axes[0,1].set_ylabel('Q')

p5 = axes[1,1].plot(pcaIQ['filtereddata_restored'][:,1],pcaIQ['filtereddata_restored'][:,12],'.')
axes[1,1].axis('equal')
#axes[1,1].set_xlabel('I')
#axes[1,1].set_ylabel('Q')

p6 = axes[2,1].plot(pcaIQ['filtereddata_restored'][:,2],pcaIQ['filtereddata_restored'][:,13],'.')
axes[2,1].axis('equal')
#axes[2,1].set_xlabel('I')
#axes[2,1].set_ylabel('Q')

#fig.subtitle('PCA: IQ Res')

#handles, labels = axes.get_legend_handles_labels()
#fig.legend(handles,labels,bbox_to_anchor=(.95,.6),borderaxespad=0.)
fig.tight_layout()

#%%
"""
#----------------------------------------------------------------------------------------------
# This section just makes test plots.
#----------------------------------------------------------------------------------------------

plt.figure('FULL');plt.plot(signal.fftconvolve(pcaIQ['initialdata'][:,2],pcaIQ['initialdata'][:,13], 'full'),'.',label='Res 3');plt.plot(signal.fftconvolve(pcaIQ['initialdata'][:,1],pcaIQ['initialdata'][:,12], 'full'),'.',label='Res 2');plt.plot(signal.fftconvolve(pcaIQ['initialdata'][:,0],pcaIQ['initialdata'][:,11], 'full'),'.',label='Res 1');plt.plot(np.max(signal.fftconvolve(pcaIQ['initialdata'][:,2],pcaIQ['initialdata'][:,13],'full')),np.where(np.max(signal.fftconvolve(pcaIQ['initialdata'][:,2],pcaIQ['initialdata'][:,13],'full')))[0],'*');plt.legend()
plt.figure('FULL Blind');plt.plot(signal.fftconvolve(pcaIQblind['initialdata'][:,2],pcaIQblind['initialdata'][:,18], 'full'),'.',label='Bin 3');plt.plot(signal.fftconvolve(pcaIQblind['initialdata'][:,1],pcaIQblind['initialdata'][:,17], 'full'),'.',label='Bin 2');plt.plot(signal.fftconvolve(pcaIQblind['initialdata'][:,0],pcaIQblind['initialdata'][:,16], 'full'),'.',label='Bin 1');plt.plot(np.max(signal.fftconvolve(pcaIQblind['initialdata'][:,2],pcaIQblind['initialdata'][:,13],'full')),np.where(np.max(signal.fftconvolve(pcaIQblind['initialdata'][:,2],pcaIQblind['initialdata'][:,13],'full')))[0],'*');plt.legend()

plt.figure('Eigenvector of PCA Component')
plt.plot(pcaIQblind['evecs'][:,0],'.-',label='Bin 1')
plt.plot(pcaIQblind['evecs'][:,1],'.-',label='Bin 2')
plt.plot(pcaIQblind['evecs'][:,2],'.-',label='Bin 3')
plt.legend()


test= np.ones((100,))/100
plt.figure(3);plt.plot(pcaIQblind['rotateddata'][:,0])
plt.figure(3);plt.plot(np.convolve(pcaIQblind['rotateddata'][:,0],test,'same'))


#pathtodir = '/Users/tashaleebillings/Desktop/'

plt.figure('Covariance Plot')
plt.clf()
plt.subplot(311)
plt.title('I')
plt.imshow(pcaI['covmatrix'],aspect='auto',interpolation='None')
plt.subplot(312)
plt.title('Q')
plt.imshow(pcaQ['covmatrix'],aspect='auto',interpolation='None')
plt.subplot(313)
plt.title('IQ')
plt.imshow(pcaIQ['covmatrix'],aspect='auto',interpolation='None')
plt.colorbar()

plt.tight_layout()
#plt.savfig(pathtodir+'CovPlot.png')
#plt.show()

plt.figure('Correlation Plot')
plt.clf()
plt.subplot(311)
plt.title('I')
plt.imshow(pcaI['corrmatrix'],aspect='auto',interpolation='None')
plt.subplot(312)
plt.title('Q')
plt.imshow(pcaQ['corrmatrix'],aspect='auto',interpolation='None')
plt.subplot(313)
plt.title('IQ')
plt.imshow(pcaIQ['corrmatrix'],aspect='auto',interpolation='None')
plt.colorbar()

plt.tight_layout()
#plt.savfig(pathtodir+'CorPlot.png')
plt.show()

#----------------------------------------------------------------------------------------------
# Made up Data Section. I used this section to test my code against Steve's IDL
#----------------------------------------------------------------------------------------------

# Generated Dataset 1
N = 10000
t = np.linspace(0,1,num=N)
x1 = np.sin(2.*np.pi*t) #np.random.randn(N)
x1 -= x1.mean()
x2 = np.cos(2.*np.pi*t) #np.random.randn(N)
x2 -= x2.mean()
c = np.random.randn(N)
a = 10
com = a*c
x1 += com
x2 += com

# Initial Data
#A = np.array([[x1],[x2]]).squeeze()
A = np.vstack((x1,x2))


# Generated Dataset 2
nvec = 1000
gx = 10.    
gy = 10.   
var_x = 1. 
var_y = 1. 
c = np.random.randn(nvec)
x = gx*c + np.sqrt(var_x)*np.random.randn(nvec)
y = gy*c + np.sqrt(var_y)*np.random.randn(nvec)
dataAdjust = (np.vstack((x,y))).transpose()

# Generated Dataset 3
np.random.seed(1)
x = np.linspace(0, 1, 500)
h = 0.01
C = np.exp(-0.5 * (x - x[:, None]) ** 2 / h ** 2)
y = 0.8 + 0.3 * np.random.multivariate_normal(np.zeros(len(x)), C)
w = np.zeros_like(x)
w[(x > 0.12) & (x < 0.28)] = 1
y_norm = np.convolve(np.ones_like(y), w, mode='full')
valid_indices = (y_norm != 0)
y_norm = y_norm[valid_indices]
y_w = np.convolve(y, w, mode='full')[valid_indices] / y_norm
x_w = np.convolve(x, w, mode='full')[valid_indices] / y_norm
y_fft = np.fft.fft(y)
w_fft = np.fft.fft(w)
yw_fft = y_fft * w_fft
yw_final = np.fft.ifft(yw_fft)

pca = pca_analysis(data=dataAdjust, normalize_stdev=1,absolute=0,nrem=0,tol=1.e-12)
print('evals =', pca['evals'])
print('evecs =', pca['evecs'])
print('PCA Filter = ', pca['filtereddata'])
print('Initial Data = ', pca['initialdata'])
print(pca['filtereddata'].shape)
print(pca['initialdata'].shape)
print(pca['covmatrix'].shape)

"""















