# C. Hales, NRAO
# 5 Feb 2020
# python3
#
# calculate linearly polarized differential counts from
# https://ui.adsabs.harvard.edu/abs/2014MNRAS.440.3113H
# see in particular II pages 3124 and 3127-28
#
# convolve H03 total intensity counts with fracpol distribution
#
# if this code helps you, please cite the paper above

import numpy as np
import matplotlib.pyplot as plt

def H03(x):
    lgx = np.log10(x)
    a0  =  0.859
    a1  =  0.508
    a2  =  0.376
    a3  = -0.049
    a4  = -0.121
    a5  =  0.057
    a6  = -0.008
    # dN/dI (not euclidean normalized)
    return 10**(a0 + a1*lgx + a2*lgx**2 + a3*lgx**3 + a4*lgx**4 +\
                a5*lgx**5 + a6*lgx**6)/(x/1e3)**2.5

lenx = 500
x    = np.logspace(np.log10(0.01),np.log10(1100),lenx)  # units: mJy
hop  = H03(x)
#plt.loglog(x,hop)
#plt.grid()

# get polarization counts curve by convolving with fracpol distribution
mStep = 0.001
m     = np.arange(0.001,1+mStep,mStep)
mcent = np.log10(0.04)
mwdth =  0.3
mamp  = 1/(np.sqrt(2*np.pi)*m*np.log(10)*mwdth)
mdist = mamp*np.exp(-(np.log10(m)-mcent)**2/(2*mwdth**2))
#plt.semilogx(m,mdist)
#plt.grid()

hopL = hop * 0
for i in range(lenx):
    print (str(i+1)+'/'+str(lenx))
    for j in range(len(m)):
        hopL[i] = hopL[i] + H03(x[i]/m[j])/m[j]*mdist[j]*mStep

# plot euclidean normalized differential counts
plt.loglog(x,hop*(x/1e3)**(2.5),x,hopL*(x/1e3)**(2.5))
plt.grid()
