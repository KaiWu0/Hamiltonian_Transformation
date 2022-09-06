#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import sys
#sys.path.append(os.getenv('HOME')+'/script')
import cal_band

def plotband(path,ax,color='r'):
    bandinfo=cal_band.band_data(path)
    k=bandinfo['k']; eigenvalues=bandinfo['eigenvalues']; ktick=bandinfo['ktick'];
    klabel=bandinfo['klabel']; fermi=bandinfo['fermi']; software=bandinfo['software']
    fermi=eigenvalues[0,3]
    eigenvalues-=fermi
    l1=plt.plot(k,eigenvalues,color,linewidth=1)
    plt.xticks(ktick,klabel)
    for i in ktick[:-1]: plt.axvline(x=i,color='k',linewidth=0.5)
    plt.xlim(ktick[0],ktick[-1])
    #plt.ylim(-2.5,2.5)
    plt.ylabel('$E-E_f$ (eV)')
    if software in ('vasp','qe'):
        #plt.axhline(y=0,color='k',linewidth=0.5)
        pass
    elif software=='phon':
        plt.axhline(y=0,color='k',linewidth=0.5,ls='--')
        #plt.ylim(np.min(),np.max(eigenvalues)+0.5)
        plt.ylabel('Frequency (THz)')
    return eigenvalues

if __name__ == "__main__":
    fig,ax=plt.subplots(figsize=(5,3.1))
    #plt.figure(figsize=(15,9))
    e0=plotband('ht',ax,'r')
    e1=plotband('wannier',ax,'b')
    e2=plotband('nscf',ax,'g')
    #plt.savefig('band.png', bbox_inches='tight')
    plt.show()
