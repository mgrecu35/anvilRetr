import netCDF4 as nc
import numpy as np
import lidarSim as lidSim
from kazrRet import *

import matplotlib.pyplot as plt
with nc.Dataset('active_obs_2020-08-11_05:10:00_subset.nc','r') as f:
    zKu=f.variables['zKu'][:]
    zKa=f.variables['zKa'][:]
    dn=f.variables['dn'][:]
    qv=f.variables['qv'][:]
    height=f.variables['height'][:]
    press=f.variables['press'][:]
    temp=f.variables['temp'][:]
    tempf=f.variables['tempf'][:]

dn_mean=np.zeros((70),float)
dn_count=np.zeros((70),float)
for i,zKu1 in enumerate(zKu):
    a=np.nonzero(zKu1>-20)
    for k in a[0]:
        dn_mean[k]+=dn[i,k]
        dn_count[k]+=1

a=np.nonzero(dn_count>10)
dn_mean[a]/=dn_count[a]
plt.plot(dn_mean[a],a[0])

dn0=np.zeros((70),float)
dn0[0:10]=-3
dn0[10:41]=-3+5*np.arange(31)/30
dn0[40:]=2

#plt.plot(dn0,range(70))
ireturn=1
freqs=['94','186.31', '325.15','660.00']
fisot=2.7
lidar=lidSim
umu=np.cos(53/180*np.pi)
iwcL=[]
tbL=[]
zKuL=[]
import tqdm
for i in tqdm.tqdm(range(zKu.shape[0])):
    zKu1=zKu[i]
    a=np.nonzero(zKu1>-20)
    iwc1=np.zeros((70),float)
    kextI=np.zeros((70,4),float)
    kscaI=np.zeros((70,4),float)
    asymI=np.zeros((70,4),float)
    for k in a[0]:
        dn01=dn0[k]
        if zKu1[k]-10*dn01<=zST[0,0]:
            dn01=(zKu1[k]-zST[0,0])/10
        if zKu1[k]-10*dn01>=zST[0,-1]:
            dn01=(zKu1[k]-zST[0,-1])/10
        ifind = lidSim.bisection2(zST[0,:],zKu1[k]-10*dn01)
        if ifind>=400:
            ifind=399
        iwc_ret=iwcST[0,ifind]*10**dn01
        iwc1[k]=iwc_ret
        kextI[k,:]=kextST[-4:,ifind]*10**dn01
        kscaI[k,:]=kscaST[-4:,ifind]*10**dn01
        asymI[k,:]=gST[-4:,ifind]
    nz=70
    
    kexttot_atm=np.zeros((nz,4),float)
    kexttot_atm=np.zeros((nz,4),float)
    emis=0.8+np.random.random()*0.1
    ebar=0.8+np.random.random()*0.1
    rho1d=press[i]/(287.*temp[i])
    for ik,freq in enumerate(freqs):
        for k in range(nz):
            absair,abswv = lidar.gasabsr98(float(freq),temp[i][k],rho1d*qv[i][k]*1e-3,press[i][k],ireturn)
            kexttot_atm[k,ik]=(absair+abswv)
    kext=kexttot_atm[:,:]+kextI[:,:]
    salb=kscaI[:,:]
    asym_fact=kscaI[:,:]*asymI
    a=np.nonzero(salb>1e-4)
    asym_fact[a]/=salb[a]
    salb/=kext
    ik=0
    tb1=[]
    for ik in range(4):
        tbout = lidar.radtran(umu,temp[i][0],temp[i],height[i],kext[:-1,ik],salb[:-1,ik],asym_fact[:-1,ik],fisot,emis,ebar)
        tb1.append(tbout)
    iwcL.append(iwc1)
    zKuL.append(zKu1)
    tbL.append(tb1)
    #break

import xarray as xr
ds=xr.Dataset({'tb':(['nt','n4'],np.array(tbL)),\
               'zKu':(['nt','nz'],np.array(zKuL)),
               'iwc':(['nt','nz'],np.array(iwcL))})
ds.to_netcdf('radiometer_retrieved_2020-08-11_05:10:00_subset.nc',\
             encoding={'tb':{'zlib':True},'zKu':{'zlib':True},'iwc':{'zlib':True}})

