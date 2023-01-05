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

nt=temp.shape[0]

tempL=[]
qvL=[]
rhoL=[]
pressL=[]
for i in range(nt):
    a=np.nonzero(temp[i]>273.15)
    hg=height[0,a[0][-1]]-14*0.240+np.arange(64)*0.240
    hgf=height[0,a[0][-1]]-14*0.240-0.120+np.arange(65)*0.240
    tempg=np.interp(hg,height[i,:],temp[i,:])
    tempgf=np.interp(hgf,height[i,:],temp[i,:])
    pressg=np.interp(hg,height[i,:],press[i,:])
    qvg=np.interp(hg,height[i,:],qv[i,:])
    rho=press[i]/287/temp[i]
    rhog=np.interp(hg,height[i,:],rho)
    tempL.append(tempgf)
    pressL.append(pressg)
    rhoL.append(rhog)
    qvL.append(qvg)
    #break

#convert tempL, rhoL, qvL and pressL to numpy arrays
tempL=np.array(tempL)
rhoL=np.array(rhoL)
qvL=np.array(qvL)
pressL=np.array(pressL)
#concatenate tempL and qvL along the second axis
tempqv=np.concatenate((tempL,qvL),axis=1)
#normalize tempqv using sklearn standardScaler
import sklearn  

from sklearn.preprocessing import StandardScaler
# define standard scaler
scaler = StandardScaler()

tempqv=scaler.fit_transform(tempqv)
#tempqv=(tempqv-np.mean(tempqv,axis=0))/np.std(tempqv,axis=0)
#cluster tempqv into 30 clusters using kmeans
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=30, random_state=0).fit(tempqv)
tempCL=[]
qvCL=[]
rhoCL=[]
pressCL=[]

for i in range(30):
    a=np.nonzero(kmeans.labels_==i)
    tempCL.append(np.mean(tempL[a],axis=0))
    qvCL.append(np.mean(qvL[a],axis=0))
    rhoCL.append(np.mean(rhoL[a],axis=0))
    pressCL.append(np.mean(pressL[a],axis=0))

#plot tempCL, qvCL, rhoCL and pressCL in separate panels
fig,ax=plt.subplots(4,1,figsize=(10,10))
for i in range(30):
    ax[0].semilogx(tempCL[i][:-1],hg)
    ax[1].semilogx(qvCL[i],hg)
    ax[2].semilogx(rhoCL[i],hg)
    ax[3].semilogx(pressCL[i],hgf)

plt.show()

#save tempCL, qvCL, rhoCL and pressCL to compressed netcdf file using xarray
import xarray as xr
ds=xr.Dataset({'temp':(['nc','nz1'],tempCL),'qv':(['nc','nz'],qvCL),'rho':(['nc','nz'],rhoCL),'press':(['nc','nz'],pressCL)})
ds.to_netcdf('profiles.nc',encoding={'temp':{'zlib':True},'qv':{'zlib':True},'rho':{'zlib':True},'press':{'zlib':True}})

