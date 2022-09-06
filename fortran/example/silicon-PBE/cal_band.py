#!/usr/bin/env python
#A general program calculates band structure for all software.
#Length: angstrom. Energy: eV.
import numpy as np
import importlib
import matplotlib.pyplot as plt
import os

hartree=27.211386245988
bohr=0.529177210903

def judge_software(path):
    s=os.listdir(path)
    software=''
    if {'band.conf','band.yaml'}.issubset(s):
        software='phon'
    elif {'OUTCAR','EIGENVAL'}.issubset(s):
        software='vasp'
    elif {'scf.in','band.in'}.issubset(s):
        software='qe'
    elif {'KSSOLV_band.hdf5'}.issubset(s):
        software='kssolv'
    for f in s:
        if f[-4:]=='.win':
            software='wannier'
    assert software, 'Cannot judge which software!'
    libband=importlib.import_module('load_'+software)
    return software,libband

def band_data(path):
    #cell:   unit cell in real space.
    #kpath:  nodes of kpoint pathes, in crystal coordinate
    #klen:   length of each path
    #klabel: label of each node
    software,libband=judge_software(path)
    if software=='vasp':
        cry=libband.read_xml(path)
        cell=cry['cell']
        fermi=cry['fermi']
        eigenvalues=cry['eigenvalues']
        kinfo=libband.read_KPOINTS(path)
        kpath=kinfo['kpath']
        klen=kinfo['klen']
        klabel=kinfo['klabel']
        idx=np.cumsum(klen[:-2])-1
        for i in reversed(idx):
            eigenvalues=np.delete(eigenvalues,i,0)
        klen-=1
    elif software=='qe':
        prefix,outdir=libband.read_scf(path)
        cry=libband.read_xml(path,prefix,outdir)
        kinfo=libband.read_band(path)
        cell=cry['cell']*bohr
        fermi=cry['fermi']*hartree
        eigenvalues=cry['eigenvalues']
        kpath=kinfo['kpath']
        klen=kinfo['klen']
        klabel=kinfo['klabel']
        if os.path.exists(path+'/band.txt'):
            print('Read Hamiltonian transformation eigenvalues!')
            eigenvalues=np.loadtxt(path+'/band.txt').transpose()
    elif software=='wannier':
        seedname=libband.get_seedname(path)
        nk,kpoints=libband.read_kpt(path,seedname)
        ktick,klabel=libband.read_labelinfo(path,seedname)
        k,eigenvalues=libband.read_eigenvalues(path,nk,seedname)
        bandinfo={'k':k,'eigenvalues':eigenvalues,'ktick':ktick,'klabel':klabel,'fermi':0,'software':software}
        return bandinfo
    else:
        print('Software not supported yet.')
    rec=2*np.pi*np.asarray(np.mat(cell).I.T)
    kpath=np.dot(kpath,rec)
    klen1=klen.copy()
    klen1[klen1<1]=1
    kpos=np.zeros([np.sum(klen1[:-1])+1,3])
    k=np.zeros([np.sum(klen1[:-1])+1])
    kcumsum=np.insert(np.cumsum(klen1[:-1]),0,0)
    label=[]
    ktick=[]
    idx=-1
    for i in range(np.shape(klen)[0]-1):
        if klen[i]>0:
            label.append(klabel[i])
            ktick.append(k[idx+1])
        else:
            idx=kcumsum[i]
            kpos[idx,:]=kpath[i,:]
            k[idx+1]=k[idx]
            label[-1]=label[-1]+','+klabel[i]
        for j in range(klen[i]):
            idx=kcumsum[i]+j
            kpos[idx,:]=kpath[i,:]+(kpath[i+1,:]-kpath[i,:])*j/klen[i]
            k[idx+1]=k[idx]+np.linalg.norm(kpath[i+1,:]-kpath[i,:])/klen[i]
    kpos[-1,:]=kpath[-1,:]
    label.append(klabel[-1])
    ktick.append(k[-1])
    bandinfo={'k':k,'eigenvalues':eigenvalues,'kpos':kpos,'ktick':ktick,'klabel':label,'fermi':fermi,'software':software}
    return bandinfo
