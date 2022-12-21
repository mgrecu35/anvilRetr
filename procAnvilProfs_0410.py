from kazrRet import *
from netCDF4 import Dataset
import numpy as np
from scipy.special import gamma as gam
import numpy as np
import lidarSim as lidar

def nw_lambd(swc,nc,mu):
    rhow=1e6
    lambd=(nc*rhow*np.pi*gam(4+mu)/gam(1+mu)/6.0/swc)**(0.333)  # m-1
    n0=nc*lambd/gam(1+mu) # m-4
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    return n0,lambd
ic=-1
from scipy.ndimage import gaussian_filter
import lidarSim as lidSim
dmCoeffs=np.polyfit(zKuG,np.log(dmG),1)
zCoeffs=np.polyfit(np.log10(gwc),zKuG,1)
import matplotlib
import xarray as xr

fh=Dataset("anvProfs_2020-08-11_04:10:00_subset.nc")
swc=fh.variables['swc'][:]
gwc=fh.variables['gwc'][:]
iwc=fh.variables['iwc'][:]
ns=fh.variables['ns'][:]
ng=fh.variables['ng'][:]
ni=fh.variables['ni'][:]
z=fh.variables['z'][:]
t=fh.variables['t'][:]
prs=fh.variables['prs'][:]
qv=fh.variables['qv'][:]
rho=fh.variables['rho'][:]
fh.close()
mu=2

def calcZ(swc,ns,mu,i,zS1):
    a=np.nonzero(swc[i]>1e-3)
    ns1,lambd1=nw_lambd(swc[i][a],ns[i][a],mu)
    dn=np.log10(ns1*1e3/0.08e8)
    dn[dn<-3]=-3
    dn[dn>3]=3
    dn+=0
    #print(a[0])
    #print(swc[i][a])
    #print(dn)
    kext=np.zeros((55,4),float)
    ksca=np.zeros((55,4),float)
    asym=np.zeros((55,4),float)
    dm=np.zeros((55),float)
    for k,swc1 in enumerate(swc[i][a]):
        ifind = lidSim.bisection2(iwcST[0,:],swc1/10**dn[k])
        if ifind>=400:
            ifind=399
            dn[k]=np.log10(swc1/iwcST[0,399])
        zs1=zST[0,ifind]+10*dn[k]
        zS1[a[0][k]]=zs1
        kext[a[0][k],:]=kextST[-4:,ifind]*10**dn[k]
        ksca[a[0][k],:]=kscaST[-4:,ifind]*10**dn[k]
        asym[a[0][k],:]=gST[-4:,ifind]
        dm[a[0][k]]=dmST[0,ifind]
    return kext,ksca,asym,dm

zKuL=[]
nt=0
iwp=(swc+gwc+iwc)*(z[:,1:]-z[:,:-1])
iwp=iwp.sum(axis=-1)
iwpL=[]
npixL=[]
aL=[]
npart=4
nrefl=2
ice_type=1
undef=0.0
pnormL=[]
iwcL=[]
zKuL=[]

def get_pnorm():
    nz3=60
    dr=0.25
    h1=z[i][0]+dr/2+np.arange(nz3)*dr
    h1f=z[i][0]+np.arange(nz3+1)*dr
    zm1=0.5*(z[i][1:]+z[i][:-1])
    temp1=np.interp(h1,zm1,t[i])
    pres1=np.interp(h1,zm1,prs[i])
    presf1=np.interp(h1f,zm1,prs[i])
    rho1=np.interp(h1,zm1,rho[i])
    iwc1=np.interp(h1,zm1,swc[i]+gwc[i]+iwc[i])
    dm_tot=(swc[i]*dmS+gwc[i]*dmG+iwc[i]*dmI)/(swc[i]+gwc[i]+iwc[i]+1e-4)
    dm_ice=np.interp(h1,zm1,dm_tot)
    q_lsice1=iwc1/rho1*1e-3
    #iwc_lidar[:,i]=iwc2d[i,::3].data
    q_lsice=q_lsice1[np.newaxis,:]
    pres1=pres1[np.newaxis,:]
    presf1=presf1[np.newaxis,:]
    temp1=temp1[np.newaxis,:]
    q_lsliq=np.zeros((1,nz3),float)
    ls_radice=dm_ice/2*1e-3
    ls_radice=ls_radice[np.newaxis,:]
    ls_radliq=np.zeros((1,nz3),float)
    q_cvice=np.zeros((1,nz3),float)
    cv_radice=np.zeros((1,nz3),float)
    q_cvliq=np.zeros((1,nz3),float)
    cv_radliq=np.zeros((1,nz3),float)
    temp=temp1
    pres=pres1
    presf=presf1
    zKuL.append(np.interp(h1,zm1,zKu))
    pmol,pnorm,pnorm_perp_tot,\
        tautot,betatot_liq,\
        betatot_ice,\
        betatot,refl, \
        zheight,\
        beta_mol, tau_mol,\
        alpha= lidar.lidar_simulator(npart,nrefl,undef,\
                                     pres1,presf1,\
                                     temp1,
                                     q_lsliq,q_lsice,\
                                     q_cvliq,\
                                     q_cvice,\
                                     ls_radliq,\
                                     ls_radice,\
                                     cv_radliq,cv_radice,\
                                     ice_type)
    pnormL.append(pnorm)
    iwcL.append(iwc1)
    

freqs=['94','186.31', '325.15','660.00']
ireturn=1
fisot=2.7
emis=0.9
ebar=0.9
umu=np.cos(53/180*np.pi)
tbL=[]
import tqdm
for i in tqdm.tqdm(range(qv.shape[0])):
    zS1=np.zeros((55),float)-99
    zG1=np.zeros((55),float)-99
    zI1=np.zeros((55),float)-99
    kextS,kscaS,asymS,dmS=calcZ(swc,ns,mu,i,zS1)
    kextG,kscaG,asymG,dmG=calcZ(gwc,ng,mu,i,zG1)
    kextI,kscaI,asymI,dmI=calcZ(iwc,ni,mu,i,zI1)
    zKu=np.log10(10.**(0.1*zS1)+10.**(0.1*zI1)+10.**(0.1*zG1))*10
    a1=np.nonzero(zKu>8)
    if len(a1[0]>2):
        iwpL.append(iwp[i])
        npixL.append(len(a1[0]))
        nt+=1
        aL.append(i)
    #continue
    nz=55
    kexttot_atm=np.zeros((nz,4),float)
    kexttot_atm=np.zeros((nz,4),float)
    for ik,freq in enumerate(freqs):
        for k in range(nz):
            absair,abswv = lidar.gasabsr98(float(freq),t[i][k],rho[i][k]*qv[i][k],prs[i][k],ireturn)
            kexttot_atm[k,ik]=(absair+abswv)
    kext=kexttot_atm[:,:]+kextS[:,:]+kextG[:,:]+kextI[:,:]
    salb=kscaS[:,:]+kscaG[:,:]+kscaI[:,:]
    asym_fact=kscaS[:,:]*asymS+kscaG[:,:]*asymG+kscaI[:,:]*asymI
    a=np.nonzero(salb>1e-4)
    asym_fact[a]/=salb[a]
    salb/=kext
    ik=0
    tb1=[]
    for ik in range(4):
        tbout = lidar.radtran(umu,t[i][0],t[i],z[i][:-1],kext[:-1,ik],salb[:-1,ik],asym_fact[:-1,ik],fisot,emis,ebar)
        tb1.append(tbout)
    tbL.append(tb1)
    get_pnorm()
    
        
aL=np.array(aL)
ds=xr.Dataset({'tb':(['nt','n4'],np.array(tbL))})
ds.to_netcdf('radiometer_obs_2020-08-11_04:10:00_subset.nc',encoding={'tb':{'zlib':True}})

ds=xr.Dataset({'iwc':(['nt','nz'],np.array(iwcL)),\
               'pnorm':(['nt','nz'],np.array(pnormL)[:,0,:]),\
               'zKu':(['nt','nz'],np.array(zKuL))})
ds.to_netcdf('obs_2020-08-11_04:10:00_subset.nc',encoding={'iwc':{'zlib':True},'pnorm':{'zlib':True},'zKu':{'zlib':True}})

def write_subset():
    pass
    ds=xr.Dataset({'swc':(['nt','nz'],swc[aL]),\
                   'gwc':(['nt','nz'],gwc[aL]),\
                   'iwc':(['nt','nz'],iwc[aL]),\
                   'ns':(['nt','nz'],ns[aL]),\
                   'ng':(['nt','nz'],ng[aL]),\
                   'ni':(['nt','nz'],ni[aL]),\
                   'z':(['nt','nz1'],z[aL]),\
               't':(['nt','nz'],t[aL]),\
                   'prs':(['nt','nz'],prs[aL]),\
                   'qv':(['nt','nz'],qv[aL]),\
                   'rho':(['nt','nz'],rho[aL])})
    
    
    ds.to_netcdf('anvProfs_2020-08-11_04:10:00_subset.nc',encoding={'swc':{'zlib':True},'gwc':{'zlib':True},'iwc':{'zlib':True},'ns':{'zlib':True},'ng':{'zlib':True},'ni':{'zlib':True},'z':{'zlib':True},'t':{'zlib':True},'prs':{'zlib':True},'qv':{'zlib':True},'rho':{'zlib':True}})


#iwc2d[ik,k]=swc[ifind]*10**dn1
#zKu=zKuS[ifind]
#zKu_2d[ik,k]=zKu+10*dn1
#dm_ice[ik,k]=np.exp(dmCoeffs[0]*(zKu)+dmCoeffs[1])
#ibin2=lidSim.bisection2(dmST[-1,:],dm_ice[ik,k])
#iwc2=iwcST[-1,ibin2]*10**dn1
#kexttot_2d[ik,k,:]=kextST[-4:,ibin2]*10**dn1
#kscatot_2d[ik,k,:]=kscaST[-4:,ibin2]*10**dn1
#asymtot_2d[ik,k,:]=gST[-4:,ibin2]



    
