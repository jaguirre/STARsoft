
# coding: utf-8

# In[1]:

import numpy as np
import h5py

#%%
def unpack(mdict,hdf5_file,grp):
    for key in mdict.keys():
        print(key)
        if grp == 0: g0 = hdf5_file.create_group(str(key))
        else: g0 = grp.create_group(str(key))
        
        if type(mdict[key]) == dict: unpack(mdict[key],hdf5_file,g0)
        else: g0.create_dataset(str(key),data=mdict[key])
    return hdf5_file

def dict_to_hdf5(mdict,filename):
    f = h5py.File(filename,'w')
    hdf5_file = unpack(mdict,f,0)
    hdf5_file.close()
    
# In[14]:

#test_dict = {0:{'s21':np.zeros(10,dtype='complex128'),'freqs':np.arange(10)},1:{'s21':np.zeros(10,dtype='complex128'),'freqs':np.arange(10)}}
#test_dict = measurement
test_dict = {0:{'n_sweeps':3,0:{'Qi':1e5,'f0':.355},1:{'Qi':1.44e5,'f0':.3554},2:{'Qi':1.5e5,'f0':.356}},
             1:{'n_sweeps':3,0:{'Qi':2e5,'f0':.365},1:{'Qi':2.44e5,'f0':.3654},2:{'Qi':2.5e5,'f0':.366}},
             2:{'n_sweeps':4,0:{'Qi':3e5,'f0':.375},1:{'Qi':3.44e5,'f0':.3754},2:{'Qi':3.5e5,'f0':.376}}}

filename = '20180806_testfile.hdf5'
dict_to_hdf5(test_dict,filename)

test_redict = {}

#%%
def repack(hdf5_file):
    mdict = {}
    
    def scan_vals(g, m):
        for k, v in g.items():
            if isinstance(v, h5py.Group):
                m[v.name] = {}
                scan_vals(v, m[v])
            elif isinstance(v, h5py.Dataset):
                m[v.name] = v
    with h5py.File(hdf5_file, 'r') as f:
        scan_vals(f,mdict)
    
    return mdict

dictfile = 'C:/Users/Alyssa/Penn Google Drive/Penn & NSTRF/Caltech Devices/test/CD012/testresdict.hdf5'
mtestdict = testdict[0][0]
dict_to_hdf5(mtestdict,dictfile)

m = repack(dictfile)
mopen = h5py.File(dictfile)

#%%
def scan_hdf5(path, recursive=True, tab_step=2):
    def scan_node(g, tabs=0):
        print(' ' * tabs, g.name)
        for k, v in g.items():
            if isinstance(v, h5py.Dataset):
                print(' ' * tabs + ' ' * tab_step + ' -', v.name)
            elif isinstance(v, h5py.Group) and recursive:
                scan_node(v, tabs=tabs + tab_step)
    with h5py.File(path, 'r') as f:
        scan_node(f)
        
dictfile = 'C:/Users/Alyssa/Penn Google Drive/Penn & NSTRF/Caltech Devices/test/CD012/CD012_datadict_dictconv.hdf5'
f = h5py.File(filename,'r')

# In[15]:

#test_dict.keys()
#
#
## In[16]:
#
#test_dict[0].keys()
#
#
## In[37]:
#
#f = h5py.File('mytestfile3.hdf5', 'w')


# In[38]:
#for key in test_dict.keys():
#    print(str(key))
#    grp = f.create_group(str(key))
#    for k in test_dict[key].keys():
#        print(str(k))
#        grp.create_dataset(str(k),data=test_dict[key][k])

#for key in test_dict.keys():
#    print(str(key))
#    grp = f.create_group(str(key))
#    for k in test_dict[key].keys():
#        print(str(k))
#        grp2 =grp.create_group(str(k))
#        for k2 in test_dict[key][k].keys():
#            grp.create_dataset(str(k2),data=test_dict[key][k][k2])
#


# In[39]:

f.close()


# In[83]:

# Create a new dictionary for reading the data back into
new_dict = {}


# In[84]:

f = h5py.File('mytestfile3.hdf5', 'r')


# In[85]:

# How to just grab some data from the hdf5 file
testdata = f['0/s21'][:]


# In[86]:

# How to read back into a dictionary of the same type we wrote out
for grp in f:
    print(grp)
    new_dict[int(grp)] ={}
    for dset in f[grp]:
        print(dset)
        new_dict[int(grp)][dset] = f[grp][dset][:]


# In[87]:

f.close()


# In[88]:

testdata


# In[79]:

new_dict


# In[80]:

new_dict[0]['s21']


# In[81]:

new_dict[0]['freqs']

