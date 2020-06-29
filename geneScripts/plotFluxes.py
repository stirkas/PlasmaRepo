#!/usr/bin/env python

import numpy as np
import numpy.polynomial.polynomial as poly
import math
import matplotlib.pyplot as plt
from plotFlux import readFlux
from plotFlux import getAverageVal
from plotFlux import getData
from plotFlux import getStdDev

dataType  = 0 #Particle flux.
avgFlux   = []
stdDev    = []
ions      = ["H$^+$", "Custom", "Be$^{+4}$", "Ne$^{+10}$", "Ar$^{+15}$", "Mo$^{+31}$", "W$^{+40}$"]
dataFiles = ["./GoerlerImpurities/nrgsummary_proton.dat",     "./GoerlerImpurities/nrgsummary_custom.dat",
             "./GoerlerImpurities/nrgsummary_beryllium.dat",  "./GoerlerImpurities/nrgsummary_neon.dat",
             "./GoerlerImpurities/nrgsummary_argon.dat",      "./GoerlerImpurities/nrgsummary_molybdenum.dat",
             "./GoerlerImpurities/nrgsummary_tungsten.dat"]
qmRat     = [1/1, 6/8, 4/9.0121820, 10/20.17970, 15/39.948, 31/95.95, 40/183.84]
startT    = [.53, .64, .37, .16, .16, .13, .16] #% location to start at. For ignoring linear phase.
ltRat     = [6.96, 6.96, 6.96, 6.96, 6.96, 6.96, 6.96] #Ratio of R/L_T = omt to normalize flux different from GENE.
densRat   = [.001, .001, .001, .001, .001, .001, .001] #Ratio of n_imp/n_e to normalize flux to impurities unlike e- from GENE.

#Add some adiabatic info too.
adiabaticIons = ["H$^+$", "W$^{+40}$"]
adiabaticDataFiles = ["./GoerlerImpurities/nrgsummary_proton_ad.dat", "./GoerlerImpurities/nrgsummary_tungsten_ad.dat"]
qmRatAd   = [1/1, 40/183.84]
startTAd  = [0, 0]
ltRatAd   = [6.96, 6.96]
densRatAd = [.001, .001]
avgFluxAd = []
stdDevAd  = []

for i,ion in enumerate(ions):
   data = []
   readFlux(dataFiles[i], data) 
   [t, fluxData] = getData(dataType, data)
   startIndex = math.ceil(startT[i]*len(fluxData)) #Get index closest to startT fraction of total time.
   fluxData = fluxData[startIndex:]
   fluxData = np.array(fluxData) #Need to make a numpy array to perform multiplication by a scalar...
   fluxData = fluxData * ((1/ltRat[i])**2) * (1/densRat[i]) #Normalize flux for L_T and n_imp.
   avgFlux.append(getAverageVal(fluxData))
   stdDev.append(getStdDev(fluxData))

for i, ion in enumerate(adiabaticIons):
   data = []
   readFlux(adiabaticDataFiles[i], data) 
   [t, fluxDataAd] = getData(dataType, data)
   startIndex = math.ceil(startT[i]*len(fluxDataAd)) #Get index closest to startT fraction of total time.
   fluxDataAd = fluxDataAd[startIndex:]
   fluxDataAd = np.array(fluxDataAd) #Need to make a numpy array to perform multiplication by a scalar...
   fluxDataAd = fluxDataAd * ((1/ltRatAd[i])**2) * (1/densRatAd[i]) #Normalize flux for L_T and n_imp.
   avgFluxAd.append(getAverageVal(fluxDataAd))
   stdDevAd.append(getStdDev(fluxDataAd))

#Plot ion points and best fit line.
fig, ax = plt.subplots()
plt.scatter(qmRat, avgFlux)
plt.errorbar(qmRat, avgFlux, yerr=stdDev, linestyle="None")

for i, txt in enumerate(ions):
    ax.annotate(txt, (qmRat[i], avgFlux[i] - stdDev[i])) #Offset y for readibility.

#Generate line of best fit.
qmRat_new = np.linspace(qmRat[len(qmRat)-1]/1.5, 1, 100) #Note: Extend low end a little past last point.
coeffs = poly.polyfit(qmRat, avgFlux, 1)
ffit   = poly.polyval(qmRat_new, coeffs)
plt.plot(qmRat_new, ffit, color='r', label='KineticEl')

#Plot adiabatic ions and line of best fit.
plt.scatter(qmRatAd, avgFluxAd)
plt.errorbar(qmRatAd, avgFluxAd, yerr=stdDevAd, linestyle="None", color='b')

for i, txt in enumerate(adiabaticIons):
    ax.annotate(txt, (qmRatAd[i], avgFluxAd[i] - stdDevAd[i])) #Offset y for readibility.

#Generate line of best fit.
qmRat_newAd = np.linspace(qmRatAd[len(qmRatAd)-1]/1.5, 1, 100) #Note: Extend low end a little past last point.
coeffsAd = poly.polyfit(qmRatAd, avgFluxAd, 1)
ffit   = poly.polyval(qmRat_newAd, coeffsAd)
plt.plot(qmRat_newAd, ffit, color='orange', label='AdiabaticEl')

plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
textSize = 12
plt.xlabel("[q/m] / [q/m]$_{H^+}$", fontsize=textSize)
plt.ylabel("<$\\Gamma$$_{ES}$>",    fontsize=textSize)
plt.tight_layout()
plt.legend(loc='upper left')
plt.show()
#plt.savefig('./GoerlerImpurities/FluxPlot_AllTime.pdf')