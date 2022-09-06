#!/usr/bin/env python
#2021/5/4
import numpy as np
import pandas as pd
import os

def get_seedname(path):
    seedname=''
    for f in os.listdir(path):
        if f[-4:]=='.win':
            seedname=f[:-4]
    assert seedname, 'Error: cannot find seedname!'
    print(seedname)
    return seedname

def read_win(path,seedname=''):
    if not seedname: seedname=get_seedname(path)
    with open(path+'/'+seedname+'.win','r') as f:
        data=f.readlines()
    data=[x.split() for x in data]
    arg=dict()
    for i,l in enumerate(data):
        if not l: continue
        if len(l)>=3 and l[1]=='=' and (l[0][0] not in ('#','!')):
            arg[l[0]]=l[2] if len(l)==3 else l[2:]
        elif l[0].lower()=='begin':
            if l[1].lower()=='unit_cell_cart':
                if data[i+4][1].lower()!='end':
                    i=i+1
                arg['cell']=np.array(data[i+1:i+4],dtype='float')
    return arg

def read_wout(path,seedname=''):
    if not seedname: seedname=get_seedname(path)
    with open(path+'/'+seedname+'.wout','r') as f:
        data=f.readlines()
    data=[x.split() for x in data]
    mode=0
    i=0
    data1=[]; data2=[]; data3=[]
    while i<len(data):
        l=data[i]
        if not l:
            i+=1
            continue
        if l[0]=='Extraction':
            mode=1
            end=i+5
            while data[end]:
                end+=1
            data1=np.array([x[:5] for x in data[i+5:end]],dtype='float')
            i=end
        elif l[-1]=='SPRD':
            data2.append([data[i-1][0],l[1],l[3],l[5]])
        elif l[-1]=='DLTA':
            data3.append([data[i-2][0],l[2],l[4],l[6]])
        i+=1
    data2=np.array(data2,dtype='float')
    data3=np.array(data3,dtype='float')
    return mode,data1,data2,data3

def read_kpt(path,seedname=''):
    if not seedname: seedname=get_seedname(path)
    kpoints=np.genfromtxt(path+'/'+seedname+'_band.kpt',dtype='float',skip_header=1)
    nk=len(kpoints)
    return nk,kpoints

def read_labelinfo(path,seedname=''):
    if not seedname: seedname=get_seedname(path)
    data=pd.read_csv(path+'/'+seedname+'_band.labelinfo.dat',header=None,delim_whitespace=True)
    klabel=list(data.iloc[:,0])
    for i in range(len(klabel)):
        if klabel[i]=='G':
            klabel[i]=r'$\Gamma$'
    ktick=data.iloc[:,2].values
    return (ktick,klabel)

def read_eigenvalues(path,nk,seedname=''):
    if not seedname: seedname=get_seedname(path)
    data=np.genfromtxt(path+'/'+seedname+'_band.dat',dtype='float')
    k=data[:nk,0]
    nband=int(data.shape[0]/nk)
    assert nk*nband==data.shape[0], 'Error: number of kpoints error!'
    eigenvalues=data[:,1].reshape((nband,nk)).transpose()
    return (k,eigenvalues)

def read_hr(path,seedname=''):
    if not seedname: seedname=get_seedname(path)
    with open(path+'/'+seedname+'_hr.dat','r') as f:
        data=f.readlines()
    data=[x.split() for x in data]
    nk=int(data[1][0])
    nband=int(data[2][0])
    start=3
    while (len(data[start])!=7) or ('.' not in data[start][5]):
        start+=1
    data=np.array(data[start:],dtype='float')
    return (nk,nband,data)

def band_data(path):
    seedname=get_seedname(path)
    nk,kpoints=read_kpt(path,seedname)
    ktick,klabel=read_labelinfo(path,seedname)
    k,eigenvalues=read_eigenvalues(path,nk,seedname)
    return (k,eigenvalues,ktick,klabel,0,0)
