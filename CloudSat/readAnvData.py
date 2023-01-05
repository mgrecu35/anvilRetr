import glob
from netCDF4 import Dataset
from csIO import *
import numpy as np
import matplotlib.pyplot as plt
from kazrRet import *
import lidarSim as lidSim
fs=sorted(glob.glob("anvData/cs*nc"))
iwp1L=[]
iwp2L=[]
for f in fs:
    fh=Dataset(f)
    iwp=fh["iwp"][:]
    zw=fh["zw"][:]
    iwc=fh["iwc"][:]
    Dm=fh["re"][:]*3.75e-3
    for i,zw1 in enumerate(zw):
        if zw1.max()<0:
            continue
        a=np.nonzero(zw1>-25)
        iwc1d=np.zeros((40),float)
        for k,zw11 in enumerate(zw1[a]):
            ifind = lidSim.bisection2(dmST[0,:],Dm[i,a[0][k]])
            dnw=(zw11-zST[2,ifind])/10.0
            iwc1d[k]=iwcST[2,ifind]*10**dnw
        iwp1L.append(iwc1d.sum())
        iwp2L.append(iwp[i])

    break
