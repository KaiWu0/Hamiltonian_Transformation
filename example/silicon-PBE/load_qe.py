#!/usr/bin/env python
#2022/1/22 modify read_band
import numpy as np
import f90nml
import xml.etree.ElementTree as ET
import os.path
import sys
sys.path.append(os.getenv('HOME')+'/script')

hartree=27.211386245988
bohr=0.529177210903

def read_scf(path):
    nml=f90nml.read(path+'/scf.in')
    prefix = nml['control']['prefix']
    outdir = nml['control']['outdir']
    return(prefix,outdir)

def read_band(path):
    with open(path+'/band.in','r') as f:
        data=f.readlines()
    data=[x.split() for x in data]
    kpath=[]
    klabel=[]
    klen=[]
    for idx,l in enumerate(data):
        if len(l)==0: continue
        if l[0]=='K_POINTS':
            mode=l[1]
            num=int(data[idx+1][0])
            for i in data[idx+2:]:
                if (len(i)==0) or (i[0][0] in ('#','!')): continue
                assert len(i)==0 or len(i)>=3, 'K_POINTS error'
                kpath.append(i[:3])
                klen.append(int(i[3]))
                if len(i)>=5:
                    if i[4]=='G':
                        i[4]=r'$\Gamma$'
                    klabel.append(i[4])
                else:
                    klabel.append('')
            break
    assert mode in ('tpiba_b','crystal_b'), 'K_POINTS mode error'
    kpath=np.array(kpath,dtype='float')
    klen=np.array(klen,dtype='int')
    kinfo={'mode':mode,'kpath':kpath,'klen':klen,'klabel':klabel}
    return kinfo

def read_xml(path,prefix,outdir):
    cry=dict()
    root = ET.parse(path+'/'+outdir+'/'+prefix+'.xml').getroot()
    cell_obj = root.findall('./output/atomic_structure/cell/')
    cry['cell'] = np.array([i.text.split() for i in cell_obj],dtype='float')
    cry['funct'] = root.find('./output/dft/functional').text
    kpt_obj = root.findall('./output/band_structure/starting_k_points/k_point')
    cry['kpoints'] = np.array([i.text.split() for i in kpt_obj],dtype='float')
    cry['weight'] = np.array([i.attrib['weight'] for i in kpt_obj],dtype='float')
    ks_obj = root.findall('./output/band_structure/ks_energies')
    cry['eigenvalues'] = np.array([i[2].text.split() for i in ks_obj],dtype='float')*hartree
    cry['fermi'] = float(root.find('./output/band_structure/fermi_energy').text)
    cry['nbnd'] = int(root.find('./output/band_structure/nbnd').text)
    cry['alat'] = float(root.findall('./output/atomic_structure')[0].attrib['alat'])
    return cry
