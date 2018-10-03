
# coding: utf-8

import numpy as np
import h5py

#%%
''' Convert a (nested) python dictionary to an HDF5 file '''
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
    
#%%
''' Convert an hdf5 file to a (nested) python dictionary '''
def repack(h5file,path,mdict): 
    for key, item in h5file[path].items():
        if isinstance(item, h5py.Dataset):
            mdict[key] = item.value
        elif isinstance(item, h5py.Group):
            mdict[key] = {}
            mdict[key] = repack(h5file,path+key+'/',mdict[key])
    return mdict    

def hdf5_to_dict(filename):
    with h5py.File(filename, 'r') as h5file:
        return repack(h5file,'/',{})

#%%
''' Prints out the structure of the hdf5 file '''
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
 
