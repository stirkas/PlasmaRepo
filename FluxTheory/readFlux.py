import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

dirName = 'GoerlerImpurities'
R_I        = np.load(dirName + '/R_I.npz')['arr_0']
R_i        = np.load(dirName + '/R_i.npz')['arr_0']
flux       = np.load(dirName + '/flux.npz')['arr_0']
growthData = np.load(dirName + '/growthData.npz')['arr_0']
species    = np.load(dirName + '/species.npz')['arr_0']

kthetaRhoi = growthData[0]
omegaReal  = growthData[1]
gamma      = growthData[2]

fig,ax = plt.subplots()
for i in range(np.shape(flux)[0]):
   fluxData = flux[i][:]
   plt.plot(growthData[0], fluxData, label=species[i])

bigFont = 14
plt.xlabel("k$_{\\theta}$$\\rho_i$", fontsize=bigFont)
plt.ylabel("<$\\Gamma_{ES}$>", fontsize=bigFont)
plt.legend()
plt.show()