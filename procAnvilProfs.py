from kazrRet import *
from netCDF4 import Dataset
import numpy as np
from scipy.special import gamma as gam
import numpy as np
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

fh=Dataset("anvProfs_2020-08-11_05:10:00_subset.nc")
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
    for k,swc1 in enumerate(swc[i][a]):
        ifind = lidSim.bisection2(iwcST[0,:],swc1/10**dn[k])
        if ifind>=400:
            ifind=399
            dn[k]=np.log10(swc1/iwcST[0,399])
        zs1=zST[0,ifind]+10*dn[k]
        zS1[a[0][k]]=zs1
        #print(zs1)

zKuL=[]
nt=0
iwp=(swc+gwc+iwc)*(z[:,1:]-z[:,:-1])
iwp=iwp.sum(axis=-1)
iwpL=[]
npixL=[]
aL=[]
for i in range(qv.shape[0]):
    zS1=np.zeros((55),float)-99
    zG1=np.zeros((55),float)-99
    zI1=np.zeros((55),float)-99
    calcZ(swc,ns,mu,i,zS1)
    calcZ(gwc,ng,mu,i,zG1)
    calcZ(iwc,ni,mu,i,zI1)
    zKu=np.log10(10.**(0.1*zS1)+10.**(0.1*zI1)+10.**(0.1*zG1))*10
    zKuL.append(zKu)
    a1=np.nonzero(zKu>8)
    if len(a1[0]>2):
        iwpL.append(iwp[i])
        npixL.append(len(a1[0]))
        nt+=1
        aL.append(i)
        
aL=np.array(aL)
#ds=xr.Dataset({'swc':(['nt','nz'],swc[aL]),\
#               'gwc':(['nt','nz'],gwc[aL]),\
#               'iwc':(['nt','nz'],iwc[aL]),\
#               'ns':(['nt','nz'],ns[aL]),\
#               'ng':(['nt','nz'],ng[aL]),\
#               'ni':(['nt','nz'],ni[aL]),\
#               'z':(['nt','nz1'],z[aL]),\
#               't':(['nt','nz'],t[aL]),\
#               'prs':(['nt','nz'],prs[aL]),\
#               'qv':(['nt','nz'],qv[aL]),\
#               'rho':(['nt','nz'],rho[aL])})


#ds.to_netcdf('anvProfs_2020-08-11_05:10:00_subset.nc',encoding={'swc':{'zlib':True},'gwc':{'zlib':True},'iwc':{'zlib':True},'ns':{'zlib':True},'ng':{'zlib':True},'ni':{'zlib':True},'z':{'zlib':True},'t':{'zlib':True},'prs':{'zlib':True},'qv':{'zlib':True},'rho':{'zlib':True}})


#iwc2d[ik,k]=swc[ifind]*10**dn1
#zKu=zKuS[ifind]
#zKu_2d[ik,k]=zKu+10*dn1
#dm_ice[ik,k]=np.exp(dmCoeffs[0]*(zKu)+dmCoeffs[1])
#ibin2=lidSim.bisection2(dmST[-1,:],dm_ice[ik,k])
#iwc2=iwcST[-1,ibin2]*10**dn1
#kexttot_2d[ik,k,:]=kextST[-4:,ibin2]*10**dn1
#kscatot_2d[ik,k,:]=kscaST[-4:,ibin2]*10**dn1
#asymtot_2d[ik,k,:]=gST[-4:,ibin2]



    
