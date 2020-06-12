#!/usr/bin/env python

import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from plotFlux import readFlux
from plotFlux import getAverageVal
from plotFlux import getData
from plotFlux import getStdDev

dataType  = 0 #Particle flux.
avgFlux   = []
qmRat     = []
stdDev    = []
ions      = ["H$^+$", "Be$^{+4}$", "Ne$^{+10}$", "Ar$^{+15}$", "Mo$^{+31}$", "W$^{+40}$"]
dataFiles = ["./GoerlerImpurities/nrgsummary_proton.dat",     "./GoerlerImpurities/nrgsummary_beryllium.dat",
             "./GoerlerImpurities/nrgsummary_neon.dat",       "./GoerlerImpurities/nrgsummary_argon.dat",
             "./GoerlerImpurities/nrgsummary_molybdenum.dat", "./GoerlerImpurities/nrgsummary_tungsten.dat"]
qmRat     = [1/1, 4/9.0121820, 10/20.17970, 15/39.948, 31/95.95, 40/183.84]

for i,ion in enumerate(ions):
   data = []
   readFlux(dataFiles[i], data)
   #if (i != 3):
   #   data = data[len(data)//2:]
   [t, fluxData] = getData(dataType, data[len(data)//2:])
   avgFlux.append(getAverageVal(fluxData))
   stdDev.append(getStdDev(fluxData))

fig, ax = plt.subplots()
plt.scatter(qmRat, avgFlux)
plt.errorbar(qmRat, avgFlux, yerr=stdDev, linestyle="None")

for i, txt in enumerate(ions):
    ax.annotate(txt, (qmRat[i], avgFlux[i] - stdDev[i])) #Offset y for readibility.

qmRat_new = np.linspace(qmRat[len(qmRat)-1]/1.5, 1, 100)
coeffs = poly.polyfit(qmRat, avgFlux, 2)
ffit   = poly.polyval(qmRat_new, coeffs)
plt.plot(qmRat_new, ffit, color='r')

plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
textSize = 12
plt.xlabel("[q/m] / [q/m]$_{H^+}$", fontsize=textSize)
plt.ylabel("<$\\Gamma$$_{ES}$>",    fontsize=textSize)
plt.tight_layout()
plt.show()