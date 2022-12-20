from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
#f=Dataset('cm1out_000009.nc')
#dbz=f['dbz'][0,:,:,:]
#dbz=np.ma.array(dbz,mask=dbz<10)
#plt.pcolormesh(dbz[:,20,:],cmap='jet')
#plt.xlim(100,180)
#plt.ylim(0,60)
#plt.colorbar()
fname="output/wrfout_d04_2020-08-11_05:10:00"
#fname="wrfout_d02_2020-08-10_18:00:00"

def read_wrf(fname,it):
    f=Dataset(fname)
    qv=f['QVAPOR'][it,:,:,:]    # water vapor
    qr=f['QRAIN'][it,:,:,:]     # rain mixing ratio
    qs=f['QSNOW'][it,:,:,:]     # snow mixing ratio
    qc=f['QCLOUD'][it,:,:,:]    # cloud mixing ratio
    qg=f['QGRAUP'][it,:,:,:]   # graupel mixing ratio
    qi=f['QICE'][it,:,:,:]   # graupel mixing ratio
    ncr=f['QNRAIN'][it,:,:,:]     # rain mixing ratio
    nci=f['QNICE'][it,:,:,:]     # snow mixing ratio
    ncg=f['QNGRAUPEL'][it,:,:,:]  # graupel mixing ratio
    ncs=f['QNSNOW'][it,:,:,:]
    #z=f['z_coords'][:]/1000.             # height (km)
    th=f['T'][it,:,:,:]+300    # potential temperature (K)
    prs=f['P'][it,:,:,:]+f['PB'][it,:,:,:]  # pressure (Pa)
    T=th*(prs/100000)**0.286  # Temperature
    t2c=T-273.15
    #stop
    z=(f['PHB'][it,:,:,:]+f['PH'][it,:,:,:])/9.81/1000.
    xlat=f['XLAT'][0,:,:]
    xlong=f['XLONG'][0,:,:]
    R=287.058  #J*kg-1*K-1
    rho=prs/(R*T)
    #return qr,qs,qg,ncr,ncs,ncg,rho,z,T,f
    return qr,qs,qg,qi,qc,ncr,ncs,ncg,nci,qv,rho,z,T,f,prs,xlong,xlat
it=-1

#qr,qs,qg,ncr,ncs,ncg,rho,z,T,fh=read_wrf(fname,it)
swcL=[]
gwcL=[]
iwcL=[]
nsL=[]
ngL=[]
niL=[]
qvL=[]
zL=[]
rhoL=[]
tL=[]
prsL=[]
for it in range(6):
    qr,qs,qg,qi,qc,ncr,ncs,ncg,nci,qv,rho,z,T,fh,prs,xlong,xlat=read_wrf(fname,it)
    
    rwc=qr*rho*1e3
    rwp=rwc*(z[1:,:,:]-z[:-1,:,:])
    rwp2d=rwp.sum(axis=0)
    a=np.nonzero(rwp2d>0.1)
    
    swc=qs*rho*1e3
    gwc=qg*rho*1e3
    iwc=qi*rho*1e3
    ncs=ncs*rho
    ncg=ncg*rho
    nci=nci*rho
    
    for i1,i2 in zip(a[0],a[1]):
        swcL.append(swc[:,i1,i2])
        gwcL.append(gwc[:,i1,i2])
        iwcL.append(iwc[:,i1,i2])
        nsL.append(ncs[:,i1,i2])
        ngL.append(ncg[:,i1,i2])
        niL.append(nci[:,i1,i2])
        qvL.append(qv[:,i1,i2])
        rhoL.append(rho[:,i1,i2])
        zL.append(z[:,i1,i2])
        tL.append(T[:,i1,i2])
        prsL.append(prs[:,i1,i2])
        
import xarray as xr


ds=xr.Dataset({'swc':(['nt','nz'],swcL),\
               'gwc':(['nt','nz'],gwcL),\
               'iwc':(['nt','nz'],iwcL),\
               'ns':(['nt','nz'],nsL),\
               'ng':(['nt','nz'],ngL),\
               'ni':(['nt','nz'],niL),\
               'z':(['nt','nz1'],zL),\
               't':(['nt','nz'],tL),\
               'prs':(['nt','nz'],prsL),\
               'qv':(['nt','nz'],qvL),\
               'rho':(['nt','nz'],rhoL)})

ds.to_netcdf('anvProfs_2020-08-11_05:10:00.nc',encoding={'swc':{'zlib':True},'gwc':{'zlib':True},'iwc':{'zlib':True},'ns':{'zlib':True},'ng':{'zlib':True},'ni':{'zlib':True},'z':{'zlib':True},'t':{'zlib':True},'prs':{'zlib':True},'qv':{'zlib':True},'rho':{'zlib':True}})

