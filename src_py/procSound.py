#PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
import numpy as np
def getEnv():
    sndData=np.loadtxt("soundOK.txt")
    presS=sndData[:,0]
    hgtS=sndData[:,1]
    tempS=sndData[:,2]
    qvS=sndData[:,5]
    return presS,hgtS,tempS,qvS
