from kazrRet import *

import numpy as np
ic=-1
from scipy.ndimage import gaussian_filter
import lidarSim as lidSim
dmCoeffs=np.polyfit(zKuG,np.log(dmG),1)
zCoeffs=np.polyfit(np.log10(gwc),zKuG,1)
import matplotlib
import xarray as xr

#ifind = lidSim.bisection2(zKaS,zka1-10*dn1)
#iwc2d[ik,k]=swc[ifind]*10**dn1
#zKu=zKuS[ifind]
#zKu_2d[ik,k]=zKu+10*dn1
#dm_ice[ik,k]=np.exp(dmCoeffs[0]*(zKu)+dmCoeffs[1])
#ibin2=lidSim.bisection2(dmST[-1,:],dm_ice[ik,k])
#iwc2=iwcST[-1,ibin2]*10**dn1
#kexttot_2d[ik,k,:]=kextST[-4:,ibin2]*10**dn1
#kscatot_2d[ik,k,:]=kscaST[-4:,ibin2]*10**dn1
#asymtot_2d[ik,k,:]=gST[-4:,ibin2]



    
