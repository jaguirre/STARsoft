# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 10:01:48 2018

@author: ACMB

Functions for splicing and sorting KID measurement dictionary items

"""

#%%
import numpy as np
import copy

def dictwhere(d,key,value,result=None,path=None):
    ''' function to find the path to a particular key,value pair within a dictionary d 
        if the key value is a scalar, e.g. d[...]['T_BB'] = 10*u.K, the result will not include the initial key
        if the value is an item in a list/array, e.g. d[...]['CSD_bin_subtracted'] = 2e-16*np.power(u.Hz,-1), 
        the result will include the initial key as well as the index of the value in the list
        (leave result=None, path=None for the initial function call) '''
        
    if result == None:
        result = []
    if path == None:
        path = []
        
    for k, v in d.items():
        if k!=key: path.append(k)
        if isinstance(v, dict):
            dictwhere(v,key,value,result,path)
        if k == key:
            if isinstance(v,list) or (isinstance(v,np.ndarray) and not v.isscalar): # << not elegant because numpy scalars are weird
                for ind,val in enumerate(v):
                    if val==value: 
                        path.append(k)
                        path.append(ind)
                        result.append(copy.deepcopy(path))
                
            elif v==value: result.append(copy.deepcopy(path))
        if k!=key: path.pop()            
    return result     

#%%
from functools import reduce  
import operator

def get_by_path(root, items):
    return reduce(operator.getitem, items, root)    

def dictget(d,inds,key):
    ''' function to get a particular value in a dictionary d given the path of indices inds to get there
        (inverse function of dictwhere)
        use key = -1 if the full path is included in inds, otherwise use key as the final key to add to inds '''
        
    result = list()
    for path in inds:
        p = path.copy()
        if key!=-1: p.append(key)
        result.append(get_by_path(d,p))
    return result

#%%
def gen_dict_extract(key,d):
    ''' function to list all the values of a particular key within a (nested) dictionary d
        to just get a list, use list(gen_dict_extract(key,d)) '''
        
    if hasattr(d,'items'):
        for k, v in d.items():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result


#%%
def dictcut(d,conds):
    ''' function to return the paths of elements in the dictionary d with elements in common as given in conds.
        conds should be a list in the form of [[key1,var1],...,[keyn,varn]] as with dictwhere. '''
    biglist = []
    for key,value in conds:
        paths = dictwhere(d,key,value)
        biglist.append(paths)
    tuplelist = [[tuple(p) for p in biglist[ind]] for ind in np.arange(0,len(biglist))] # annoyingly, need to have tuples to use the set + intersection call below
    commonalities = list(set(tuplelist[0]).intersection(*tuplelist))
    commonalities = [list(t) for t in commonalities]
    return commonalities


#%%
def redict(d,paths):
    ''' function to plug back in and get a new dictionary with only specific paths of a bigger dictionary d
    '''
    newdict = {}
    for sweep,scan in paths:
        try:
            newdict[sweep][scan] = d[sweep][scan]
        except KeyError:
            newdict[sweep] = {}
            newdict[sweep][scan] = d[sweep][scan]
    return newdict




