#!/usr/bin/env python
import numpy as np
import sys
from scipy.special import erf
from numpy import exp,sqrt,pi

def fun(a,n,s,x):
    y=np.array(x)-s
    y=np.piecewise(y,[y>=0,y<=-a],[0,lambda y:y+a/2,lambda y:a*(exp(-n*n/4)-exp(-n*n*(a+2*y)*(a+2*y)/4/a/a))/2/n/sqrt(pi)/erf(n/2) +(a+2*y)*(erf(n/2)-erf(n*(0.5+y/a)))/4/erf(n/2)])
    return y

def dfun(a,n,s,x):
    y=np.array(x)-s
    y=np.piecewise(y,[y>=0,y<=-a],[0,1,lambda y: 0.5-erf(n*(0.5+y/a))/2/erf(n/2)])
    return y

def newton_inv(a,n,s,y):
    dxmax=a/10
    resmax=1e-14
    y[y>0]=0
    x=y+s
    res=1.
    niter=0
    while( (np.max(abs(res))>resmax) and (niter<50)):
        res=fun(a,n,s,x)-y
        df=dfun(a,n,s,x)
        dx=np.zeros_like(y)
        mask=(abs(df)!=0)
        dx[mask]=-res[mask]/df[mask]
        dx[dx> dxmax] =  dxmax
        dx[dx<-dxmax] = -dxmax
        x = x + dx
        niter = niter + 1
    if (np.max(abs(res))>resmax):
        print("Warning: Newton's method does not converge!")
    return x

assert len(sys.argv)>=3 and len(sys.argv)<=4,'Wrong number of arguments.\n\tUsage: python transform.py prefix forward/backward.'
arg=sys.argv[-1]
prefix=sys.argv[-2]

if arg=='forward':
    filename=prefix+'.eig'
    data=np.genfromtxt(filename,dtype='float')
    print('Warning: Here we truncate nbnd!')
    data=data[data[:,0]<=8,:]
    
    nbnd=int(np.max(data[:,0]))
    nk=int(np.max(data[:,1]))
    s=np.max(data[:,2])
    eig=np.reshape(data[:,2],[nk,nbnd])
    a=4*(np.max(eig[:,-1])-np.min(eig[:,-1]))
    power=0.5
    n=power
    y=fun(a,n,s,eig)
    
    data[:,2]=np.reshape(y,nbnd*nk)
    np.savetxt(filename, data, fmt=' %4d %4d % 17.12f')
    np.savetxt('ht_args',[a,n,s,nbnd])

elif arg=='backward':
    filename=prefix+'_band.dat'
    data=np.genfromtxt(filename,dtype='float')
    
    eig=data[:,1]
    m=np.loadtxt('ht_args')
    a=m[0]
    n=m[1]
    s=m[2]
    nbnd=m[3]
    nk=round(len(eig)/nbnd)
    
    x=newton_inv(a,n,s,eig)
    data[:,1]=x
    with open(filename,'w') as f:
        for i in range(len(eig)):
            if i%nk==0 and i>0:
                f.write("\n")
            f.write(f' {data[i,0]: 9.8E} {data[i,1]: 10.8E}\n')
else:
    print('Error, unrecognized input parameter!')
    exit()
