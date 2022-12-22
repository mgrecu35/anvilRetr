from netCDF4 import Dataset
from sklearn.neighbors import KNeighborsRegressor
import numpy as np
neigh = KNeighborsRegressor(n_neighbors=30,weights='distance')
from sklearn.model_selection import train_test_split

with Dataset("radiometer_obs_2020-08-11_04:10:00_subset.nc") as fh:
    yobsL=fh["tb"][:]

with Dataset("active_obs_2020-08-11_04:10:00_subset.nc") as fh:
    xL=fh["iwc"][:].sum(axis=-1)
    x2D=fh["iwc"][:,15:55]
    zku=fh["zKu"][:]
    pnorm=fh["pnorm"][:]

nt,nc=yobsL.shape
yobsL+=np.random.randn(nt,nc)*3
ind_train,ind_test,y_train, y_test = train_test_split(range(xL.shape[0]), x2D[:], test_size=0.33, random_state=42)
nt,nc=yobsL.shape
X_train=yobsL[ind_train,:]
X_test=yobsL[ind_test,:]
neigh.fit(X_train, y_train)

import numpy as np
import matplotlib.pyplot as plt

y_=neigh.predict(X_test)

ax=plt.subplot(111)
plt.scatter(y_test[:,30],y_[:,30])
#plt.xlim(0,16)
#plt.ylim(0,16)
ax.set_aspect('equal')

#plt.figure()
#plt.show()

nc=6
import matplotlib
from sklearn.cluster import KMeans

kmeans = KMeans(n_clusters=3, random_state=0).fit(X_train)
zcfad=np.zeros((22,55,nc),float)
bcfad=np.zeros((30,55,nc),float)
for i in range(nc):
    a=np.nonzero(kmeans.labels_==i)
    plt.figure()
    for i1 in a[0]:
        for k,z1 in enumerate(zku[i1,15:70]):
            i0=int(z1-8)
            if i0>=0 and i0<22:
                zcfad[i0,k]+=1
        for k,bscat1 in enumerate(pnorm[i1,15:70]):
            if bscat1<1e-8:
                continue
            i0=int((np.log10(bscat1)+6.8)*10)
            if i0>=0 and i0<30:
                bcfad[i0,k]+=1
                
    plt.subplot(121)
    plt.pcolormesh(zcfad[:,:,i].T,norm=matplotlib.colors.LogNorm(),cmap='jet')
    plt.subplot(122)
    plt.pcolormesh(bcfad[:,:,i].T,norm=matplotlib.colors.LogNorm(),cmap='jet')
