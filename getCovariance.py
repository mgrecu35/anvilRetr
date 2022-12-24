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

for i,zKu1 in enumerate(zKu):
    a=np.nonzero(zKu1>-20)
    iwc1=np.zeros((70),float)
    kext=np.zeros((70,4),float)
    ksca=np.zeros((70,4),float)
    asym=np.zeros((70,4),float)
    for k in a[0]:
        ifind = lidSim.bisection2(zST[0,k],zKu1[k]-10*dn0[k])
        iwc_ret=iwcST[0,ifind]*10**dn0[k]
        kext[k,:]=kextST[-4:,ifind]*10**dn0[k]
        ksca[k,:]=kscaST[-4:,ifind]*10**dn0[k]
        asym[k,:]=gST[-4:,ifind]
    break
