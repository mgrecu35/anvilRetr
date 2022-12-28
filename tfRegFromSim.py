from netCDF4 import Dataset
from sklearn.neighbors import KNeighborsRegressor
import numpy as np
neigh = KNeighborsRegressor(n_neighbors=30,weights='distance')
from sklearn.model_selection import train_test_split

with Dataset("radiometer_retrieved_2020-08-11_05:10:00_subset.nc") as fh:
    yobsL=fh["tb"][:,:]
    x2D=fh["iwc"][:,15:55]
    zku=fh["zKu"][:]
    
nt,nc=yobsL.shape
yobsL+=np.random.randn(nt,nc)*3
ind_train,ind_test,y_train, y_test = train_test_split(range(x2D.shape[0]), x2D[:], test_size=0.33, random_state=42)
nt,nc=yobsL.shape
X_train=yobsL[ind_train,:]
X_test=yobsL[ind_test,:]
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from sklearn.cluster import KMeans

neigh.fit(X_train, y_train)
ynn_=neigh.predict(X_test)

zku[zku<8]=0
zku_train=zku[ind_train,:]
zku_test=zku[ind_test,:]
from sklearn.preprocessing import StandardScaler

scalerIWC = StandardScaler()
iwc_train=y_train
iwc_test=y_test
scalerIWC.fit(iwc_train)                                 #scalerIWC is the IWC scaler
iwc_train_scaled=scalerIWC.transform(iwc_train)
iwc_test_scaled=scalerIWC.transform(iwc_test)

from sklearn.decomposition import PCA
pca = PCA(n_components=10)
pca.fit(iwc_train_scaled)

pca_train=pca.transform(iwc_train_scaled)
pca_test=pca.transform(iwc_test_scaled)

import tensorflow as tf
def model(n,m,k):
    inputs = tf.keras.Input(shape=(n,))
    x = tf.keras.layers.Dense(m, activation='relu')(inputs)
    #define output layer
    x = tf.keras.layers.Dropout(0.1)(x)
    x = tf.keras.layers.Dense(m, activation='relu')(x)
    x = tf.keras.layers.Dropout(0.1)(x)
    outputs = tf.keras.layers.Dense(k)(x)
    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    model.compile(optimizer='adam',loss='mae',metrics=['mae'])
    return model

xnn_train=np.concatenate((zku_train[:,15:55],X_train[:,:]),axis=-1)
xnn_test=np.concatenate((zku_test[:,15:55],X_test[:,:]),axis=-1)
scalerX=StandardScaler()
scalerX.fit(xnn_train)
xnn_train=scalerX.transform(xnn_train)                #scalerX is the observation scaler
xnn_test=scalerX.transform(xnn_test)
tfmodel=model(44,32,10)
scalerPCA=StandardScaler()                            #scalerPCA is the IWC PCA scaler
scalerPCA.fit(pca_train)
pca_train=scalerPCA.transform(pca_train)
pca_test=scalerPCA.transform(pca_test)
tfmodel.fit(xnn_train, pca_train, epochs=30, batch_size=32,\
          validation_data=(xnn_test, pca_test))

tfmodel.save('radar_radiom_tf_model.h5')
import pickle
docs="""
observations consist of 15:55 ku obs + 4 frequencies and are scaled by scalerX
the reference IWC is scaled by scalerIWC and then transformed into PCA(n=10) space
the associated model is radar_radiom_tf_model.h5
the radiometer-only nearest neighbor estimator is neigh
"""
pickle.dump({"scalerX":scalerX,"scalerIWC":scalerIWC,"scalerPCA":scalerPCA,"PCA":pca,"meta":docs,\
             "nneigh":neigh},\
            open("scalersAndPCA.pklz","wb"))

y_=tfmodel.predict(xnn_test)
y_=np.array(y_)
y_=scalerPCA.inverse_transform(y_)                   # unscale the PCAs

y_=pca.inverse_transform(y_)                         # from IWC PCA to scaled IWC
y_=scalerIWC.inverse_transform(y_)                   # unscale the IWC

yf_=y_.flatten()                                     # flatten the IWC

#ynnf_=ynn_.flatten()
iwcf=iwc_test.flatten()                              # flatten the test IWC
plt.figure()
ax1=plt.subplot(121)
plt.hist2d(yf_,iwcf,bins=np.arange(100)*0.01,\
           norm=matplotlib.colors.LogNorm(),cmap='jet')
x1=np.arange(101)*0.01
plt.plot(x1,x1)
plt.ylabel('Retrieved IWC [g/m$^3$]')
plt.xlabel('Reference IWC [g/m$^3$]')
plt.title('Synergistic')
ax1.set_aspect('equal')
cb1=plt.colorbar(orientation='horizontal')
cb1.ax.set_xlabel('Counts')

if 1==2:
    ax2=plt.subplot(122)
    plt.hist2d(ynnf_,iwcf,bins=np.arange(100)*0.01,\
                norm=matplotlib.colors.LogNorm(),cmap='jet')
    plt.plot(x1,x1)
    ax2.set_aspect('equal')
    plt.xlabel('Reference IWC [g/m$^3$]')
    plt.title('Radiometer only')
    cb2=plt.colorbar(orientation='horizontal')
    cb2.ax.set_xlabel('Counts')

plt.savefig('retrievals_synth_hist2d.png')

if 1==2:
    plt.figure(figsize=(6,8))
    hm=h[:,15:55].mean(axis=0)
    plt.plot(iwc_test.mean(axis=0),hm)
    plt.plot(yf_.reshape(4652,40).mean(axis=0)*ynnf_.mean()/yf_.mean(),hm)
    plt.plot(ynnf_.reshape(4652,40).mean(axis=0),hm)
    plt.xlabel('IWC [g/m$^3$]')
    plt.ylabel('Height [km]')
    plt.legend(['Reference','Synergistic','Radiometer-only'])
    plt.title('Mean profiles')
    plt.savefig('retrievals_profiles.png')
    
