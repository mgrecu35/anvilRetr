from netCDF4 import Dataset
from sklearn.neighbors import KNeighborsRegressor
import numpy as np
neigh = KNeighborsRegressor(n_neighbors=30,weights='distance')
from sklearn.model_selection import train_test_split

with Dataset("radiometer_obs_2020-08-11_04:10:00_subset.nc") as fh:
    yobsL=fh["tb2"][:,:]
    yobs2L=fh["tb"][:,:]
with Dataset("active_obs_2020-08-11_04:10:00_subset.nc") as fh:
    x2D=fh["iwc"][:,15:55]
    zku=fh["zKu"][:]
    
nt,nc=yobsL.shape
yobsL+=np.random.randn(nt,nc)*3
X_test=yobsL[:,:]
iwc_test=x2D
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
import pickle
zku_test=zku
import tensorflow as tf
d=pickle.load(open("dfnscalersAndPCA.pklz","rb"))
pca=d["PCA"]
scalerPCA=d["scalerPCA"]
scalerX=d["scalerX"]
neigh=d["nneigh"]
scalerIWC=d["scalerIWC"]
zku_test[zku_test<8]=0
xnn_test=np.concatenate((zku_test[:,15:55],X_test[:,:]),axis=-1)
xnn_test=scalerX.transform(xnn_test)
tfmodel=tf.keras.models.load_model('dfnradar_radiom_tf_model.h5')
ynn_=neigh.predict(X_test)

y_=tfmodel.predict(xnn_test)
y_=np.array(y_)
y_=scalerPCA.inverse_transform(y_)                   # unscale the PCAs

y_=pca.inverse_transform(y_)                         # from IWC PCA to scaled IWC
y_=scalerIWC.inverse_transform(y_)                   # unscale the IWC

yf_=y_.flatten()                                     # flatten the IWC

ynnf_=ynn_.flatten()
iwcf=iwc_test.flatten()                              # flatten the test IWC
plt.figure()
ax1=plt.subplot(121)
plt.hist2d(iwcf,yf_,bins=np.arange(100)*0.01,\
           norm=matplotlib.colors.LogNorm(),cmap='jet')
x1=np.arange(101)*0.01
plt.plot(x1,x1)
plt.ylabel('Retrieved IWC [g/m$^3$]')
plt.xlabel('Reference IWC [g/m$^3$]')
plt.title('Synergistic')
ax1.set_aspect('equal')
cb1=plt.colorbar(orientation='horizontal')
cb1.ax.set_xlabel('Counts')

if 1==1:
    ax2=plt.subplot(122)
    plt.hist2d(iwcf,ynnf_,bins=np.arange(100)*0.01,\
                norm=matplotlib.colors.LogNorm(),cmap='jet')
    plt.plot(x1,x1)
    ax2.set_aspect('equal')
    plt.xlabel('Reference IWC [g/m$^3$]')
    plt.title('Radiometer only')
    cb2=plt.colorbar(orientation='horizontal')
    cb2.ax.set_xlabel('Counts')

plt.savefig('retrievals_PCAsfromRet_noise_hist2d.png')

