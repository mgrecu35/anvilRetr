from netCDF4 import Dataset
from sklearn.neighbors import KNeighborsRegressor
import numpy as np
neigh = KNeighborsRegressor(n_neighbors=30,weights='distance')
from sklearn.model_selection import train_test_split

with Dataset("radiometer_obs_2020-08-11_04:10:00_subset.nc") as fh:
    yobsL=fh["tb"][:,1:]

with Dataset("active_obs_2020-08-11_04:10:00_subset.nc") as fh:
    xL=fh["iwc"][:].sum(axis=-1)
    x2D=fh["iwc"][:,15:55]
    zku=fh["zKu"][:]
    pnorm=fh["pnorm"][:]
    pmol=fh["pmol"][:]
    pnorm-=pmol
    h=fh["height"][:]
    
nt,nc=yobsL.shape
yobsL+=np.random.randn(nt,nc)*3
ind_train,ind_test,y_train, y_test = train_test_split(range(xL.shape[0]), x2D[:], test_size=0.33, random_state=42)
nt,nc=yobsL.shape
X_train=yobsL[ind_train,:]
X_test=yobsL[ind_test,:]
neigh.fit(X_train, y_train)

import numpy as np
import matplotlib.pyplot as plt

ynn_=neigh.predict(X_test)

#ax=plt.subplot(111)
#plt.scatter(y_test[:,30],y_[:,30])
#plt.xlim(0,16)
#plt.ylim(0,16)
#ax.set_aspect('equal')

#plt.figure()
#plt.show()

nc=30
import matplotlib
from sklearn.cluster import KMeans

kmeans = KMeans(n_clusters=nc, random_state=0).fit(X_train)
zcfad=np.zeros((42,55,nc),float)
bcfad=np.zeros((30,55,nc),float)
for i in [1]:#range(1):
    a=np.nonzero(kmeans.labels_>=0)

    for i1 in a[0][:]:
        for k,z1 in enumerate(zku[i1,15:70]):
            i0=int(z1+10)
            k0=int((h[i1,k+15]-4.0)/0.25)
            if i0>=0 and i0<42 and k0>=0 and k0<55:
                zcfad[i0,k0]+=1
        for k,bscat1 in enumerate(pnorm[i1,15:70]):
            if bscat1<1e-8:
                continue
            k0=int((h[i1,k+15]-4.0)/0.25)
            i0=int((np.log10(bscat1)+6.8)*10)
            if i0>=0 and i0<30 and k0>=0 and k0<55:
                bcfad[i0,k0]+=1

    fig=plt.figure()
    plt.subplot(121)
    plt.pcolormesh(-10+np.arange(42),4+np.arange(55)*0.25,\
                   zcfad[:,:,i].T,norm=matplotlib.colors.LogNorm(),cmap='jet')
    ax = plt.gca()  
    box = ax.get_position()
    plt.xlabel("Reflectivity [dBZ]")
    plt.ylabel("Height[km]")
    ax.set_position([box.x0, box.y0+0.2, box.width, box.height-0.2])

    cax1=fig.add_axes([box.x0, box.y0, box.width, 0.02])

    cb1=plt.colorbar(cax=cax1,orientation='horizontal')
    
    cb1.ax.set_title('Counts')
    plt.subplot(122)
    plt.pcolormesh(10**(-6.8+np.arange(30)*0.1),\
                   4+np.arange(55)*0.25,bcfad[:,:,i].T,norm=matplotlib.colors.LogNorm(),cmap='jet')
    plt.xscale('log')
    plt.xlabel("Lidar backscatter [m$^{-1}$sr$^{-1}$]")
    ax = plt.gca()  
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+0.2, box.width, box.height-0.2])

    cax2=fig.add_axes([box.x0, box.y0, box.width, 0.02])
    cb2=plt.colorbar(cax=cax2,orientation='horizontal')
    cb2.ax.set_title('Counts')
    plt.savefig('RadarAndLidarDistr.png')

    fig=plt.figure()

    plt.subplot(121)
    plt.pcolormesh(-10+np.arange(42),np.arange(55)+15,\
                   zcfad[:,:,i].T,norm=matplotlib.colors.LogNorm(),cmap='jet')
    ax = plt.gca()  
    box = ax.get_position()
    plt.xlabel("Reflectivity [dBZ]")
    plt.ylabel("Height[km]")
    ax.set_position([box.x0, box.y0+0.2, box.width, box.height-0.2])

    cax1=fig.add_axes([box.x0, box.y0, box.width, 0.02])

    cb1=plt.colorbar(cax=cax1,orientation='horizontal')
    
    cb1.ax.set_title('Counts')
    plt.subplot(122)
    plt.pcolormesh(10**(-6.8+np.arange(30)*0.1),\
                   np.arange(55)+15,bcfad[:,:,i].T,norm=matplotlib.colors.LogNorm(),cmap='jet')
    plt.xscale('log')
    plt.xlabel("Lidar backscatter [m$^{-1}$sr$^{-1}$]")
    ax = plt.gca()  
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+0.2, box.width, box.height-0.2])

    cax2=fig.add_axes([box.x0, box.y0, box.width, 0.02])
    cb2=plt.colorbar(cax=cax2,orientation='horizontal')
    cb2.ax.set_title('Counts')
    plt.savefig('RadarAndLidarDistrBin.png')


stop


nc=20
import matplotlib
from sklearn.cluster import KMeans

kmeans = KMeans(n_clusters=nc, random_state=0).fit(X_train)
zku[zku<-5]=-5
zku_train=zku[ind_train,:]
zku_test=zku[ind_test,:]
kgainL=[]
xL=[]
yL=[]
iwc_train=y_train
iwc_test=y_test
R=4
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
scaler.fit(iwc_train)
iwc_train_scaled=scaler.transform(iwc_train)
iwc_test_scaled=scaler.transform(iwc_test)

from sklearn.decomposition import PCA

pca = PCA(n_components=10)
pca.fit(iwc_train_scaled)

pca_train=pca.transform(iwc_train_scaled)
pca_test=pca.transform(iwc_test_scaled)

for i in range(nc):
    a=np.nonzero(kmeans.labels_==i)
    Y=np.concatenate((zku_train[a[0],15:55],X_train[a[0],:]),axis=-1)
    nobs=Y.shape[1]
    X=pca_train[a[0],:]
    X[X<1e-6]=1e-6
    nstate=X.shape[1]
    covXY=np.cov(Y.T,(X.T))
    covYY=covXY[0:nobs,0:nobs]
    covXY1=covXY[nobs:,0:nobs]
    kgainc=np.dot(covXY1,np.linalg.inv(covYY+R*np.eye(nobs)))
    kgainL.append(kgainc)
    xL.append(X.mean(axis=0))
    yL.append(Y.mean(axis=0))
    print(len(a[0]))

labels_=kmeans.predict(X_test)
y_=[]
yc_=[]
for i,l in enumerate(labels_):
    Y=np.concatenate((zku_test[i,15:55],X_test[i,:]),axis=-1)
    dy=Y-yL[l]
    a1=np.nonzero(zku_test[i,15:55]<8)
    dy[a1]=0
    y=((xL[l])+np.dot(kgainL[l],dy))
    y_.append(y)
    yc_.append(xL[l])

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
xnn_train=scalerX.transform(xnn_train)
xnn_test=scalerX.transform(xnn_test)
tfmodel=model(43,32,10)
scalerPCA=StandardScaler()
scalerPCA.fit(pca_train)
pca_train=scalerPCA.transform(pca_train)
pca_test=scalerPCA.transform(pca_test)
tfmodel.fit(xnn_train, pca_train, epochs=30, batch_size=32,\
          validation_data=(xnn_test, pca_test))

y_=tfmodel.predict(xnn_test)
y_=np.array(y_)
y_=scalerPCA.inverse_transform(y_)
#yc_=np.array(yc_)
y_=pca.inverse_transform(y_)
y_=scaler.inverse_transform(y_)
#yc_=pca.inverse_transform(yc_)
#yc_=scaler.inverse_transform(yc_)
yf_=y_.flatten()
#ycf_=yc_.flatten()
ynnf_=ynn_.flatten()
iwcf=iwc_test.flatten()
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

ax2=plt.subplot(122)
plt.hist2d(ynnf_,iwcf,bins=np.arange(100)*0.01,\
           norm=matplotlib.colors.LogNorm(),cmap='jet')
plt.plot(x1,x1)
ax2.set_aspect('equal')
plt.xlabel('Reference IWC [g/m$^3$]')
plt.title('Radiometer only')
cb2=plt.colorbar(orientation='horizontal')
cb2.ax.set_xlabel('Counts')
plt.savefig('retrievals_hist2d.png')
#plt.subplot(211)
#plt.hist2d(yf_,iwcf)

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
