import glob
from netCDF4 import Dataset
from csIO import *
import numpy as np
import matplotlib.pyplot as plt

def process(dbz,cpr_cmask,temp,iwc,re,height,iwp):
    a=np.nonzero(np.array(iwp)>0)
    b=np.nonzero(temp[a].max(axis=1)>280)
    a=a[0][b]
    temp0L=[]
    isfcPrecip=0
    iwp1L=[]
    iwp2L=[]
    iwcL=[]
    reL=[]
    zwL=[]
    for i in a[0:]:
        b1=np.nonzero(temp[i]<273.15)
        ik=b1[0][-1]
        nprecip1=np.nonzero(cpr_cmask[i,ik:100]>30)
        #print(ik,len(nprecip1[0]))
        temp0L.append(ik)
        if(len(nprecip1[0])>3):
            isfcPrecip+=1
            #print(cpr_cmask[i,ik:100])
            #print(nprecip1[0])
            dbz[i,:]=-99
        else:
            if ik+3-40>0:
                iwp1L.append(iwp[i])
                iwp2L.append(iwc[i,ik+3-40:ik+3].data.sum()*0.24)
                iwcL.append(iwc[i,ik+3-40:ik+3])
                reL.append(re[i,ik+3-40:ik+3])
                zwL.append(dbz[i,ik+3-40:ik+3])
    return iwp1L,iwp2L,zwL,iwcL,reL

fs=sorted(glob.glob("data/2019*hdf"))
fs_2B=sorted(glob.glob("CloudSat_2B/2019*hdf"))
import xarray as xr
for i,f in enumerate(fs[:-14]):
    iwc,dbz_sim,temp,re,height,lon,lat,iwp=read2C_ICE(f)
    fh=Dataset(fs_2B[i+14])
    dbz=fh['Radar_Reflectivity'][:]/100.
    cpr_cmask=fh['CPR_Cloud_mask'][:]
    iwp1L,iwp2L,zwL,iwcL,reL=process(dbz,cpr_cmask,temp,iwc,re,height,iwp)

    ds=xr.Dataset({'iwp':(['nt'],iwp1L),'zw':(['nt','nz'],zwL),\
                   'iwc':(['nt','nz'],iwcL),'re':(['nt','nz'],reL)})
    orbit=f.split('_')[1]
    fname_out='anvData/cs_profs_%s.nc'%orbit
    ds.to_netcdf(fname_out,encoding={'iwp':{'zlib':True},'zw':{'zlib':True},\
                                     'iwc':{'zlib':True},'re':{'zlib':True}})
    print(orbit,fs_2B[i+14].split('_')[2])
    #break



    
    
#plt.hist(iwp[a])
#n4=int(len(a)/4)
for i in range(-4):
    plt.figure(figsize=(10,8))
    plt.subplot(211)
    dbzm=np.ma.array(dbz,mask=dbz<-25)
    plt.pcolormesh(dbzm[a,:].T,cmap='jet',vmin=-20,vmax=25)
    plt.plot(temp0L)
    plt.xlim(i*n4,(i+1)*n4)
    plt.ylim(105,45)
    plt.colorbar()
    plt.subplot(212)
    plt.pcolormesh(cpr_cmask[a,:].T,cmap='jet')
    plt.xlim(i*n4,(i+1)*n4)
    plt.plot(temp0L)
    plt.ylim(105,45)
    plt.colorbar()
