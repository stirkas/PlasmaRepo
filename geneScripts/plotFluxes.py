#!/usr/bin/env python

import numpy as np
import numpy.polynomial.polynomial as poly
import math
import matplotlib.pyplot as plt
from plotFlux import readFlux
from plotFlux import getAverageVal
from plotFlux import getData
from plotFlux import getStdDev

#Data indices to use when gathering info.
names     = 0
files     = 1
qmRat     = 2
densRat   = 3
startT    = 4
omt       = 5
dataSet   = 6
dataSetColor = 7
dataFitColor = 8
#Normalization types
experimentalNorm = "ExpNorm"
impurityNorm     = "ImpNorm"
#Flux data types
particleFlux = 0
heatFlux     = 1
momentumFlux = 2

#Set data flags for script.
#Data type flag. Set just one true. For kinetic, adiabatic, or both electron data.
adiabaticData = False
kineticData   = False
bothData      = True
#Normalization flags. Set at most one true. None gives GENE norm.
impNormFlag = False  #Convert GENE norm to nicer imp. norm using L_T and n_imp.
expNormFlag = False  #Convert GENE norm to experimental norm.
fluxType  = particleFlux

#Data sets.
kineticDataSet =   [["H$^+$", "Custom$^{+6}$", "Be$^{+4}$", "Ne$^{+10}$", "Ar$^{+15}$", "Mo$^{+31}$", "W$^{+40}$"],
                    ["./GoerlerImpurities/nrgsummary_proton.dat",     "./GoerlerImpurities/nrgsummary_custom.dat",
                     "./GoerlerImpurities/nrgsummary_beryllium.dat",  "./GoerlerImpurities/nrgsummary_neon.dat",
                     "./GoerlerImpurities/nrgsummary_argon.dat",      "./GoerlerImpurities/nrgsummary_molybdenum.dat",
                     "./GoerlerImpurities/nrgsummary_tungsten.dat"],
                    [1/1, 6/8, 4/9.0121820, 10/20.17970, 15/39.948, 31/95.95, 40/183.84],
                    [.001, .001, .001, .001, .001, .001, .001],
                    [.53, .64, .37, .16, .16, .13, .16],
                    [6.96, 6.96, 6.96, 6.96, 6.96, 6.96, 6.96],
                    "kinetic",
                    "steelblue",
                    "red"]
adiabaticDataSet = [["H$^+$", "Ne$^{+10}$", "W$^{+40}$"],
                    ["./GoerlerImpurities/nrgsummary_proton_ad.dat", "./GoerlerImpurities/nrgsummary_neon_ad.dat",
                     "./GoerlerImpurities/nrgsummary_tungsten_ad.dat"],
                    [1/1, 10/20.17970, 40/183.84],
                    [.001, .001, .001],
                    [.114, .25, .118],
                    [6.96, 6.96, 6.96],
                    "adiabatic",
                    "red",
                    "blue"]


#Function for converting flux to SI units with DIII-D params. Thanks Scott!
def getExpFactor():
   # Parameters are all obtained from Callen Nuc. Fusion 2010 by Scott
   # This is a DIII-D H-mode before the onset of an ELM
   R0=1.7 # Magnetic axis, meters
   R=2.2 # location of cold impurities, near separatrix
   Bt0= 2. # Teslas, B at R0
   mu = 1. # m/mp for bulk ions
   Z = 1. # This is Z for the bulk ions
   Ti = 1000. # bulk ion temp in eV at cold impurity location
   c = 0.001 # c=nI/n0, concentration of impurity

   # Flux = Flux_GENE x n0 vt (rhoi/R)^2
   # R=L_ref, n0=n_ref, vt=vt_ref, rhoi_ref
   # vt=979.0/sqrt(mu) x sqrt(Ti)
   # I use mu=2 in this formula because tokamaks typically run with D
   # rhoi=1.02 x sqrt(mu Ti)/(Z B)
   # factor= Flux/(nI*Flux_GENE) = 1/c vt (rhoi/R)^2

   # a little tricky... We want v_pinch.  Flux = v_pinch x nI
   # Flux_GENE = v_pinch x n0.  v_pinch = Flux_GENE/(nI/n0) x vt (rhoi/R)^2

   B = Bt0*R0/R*10000. # Plasma formulary formulas use Gauss
   rhoi = 1.02*np.sqrt(mu*Ti)/(Z*B)
   vt=9790.*np.sqrt(Ti/mu)

   factor = vt/c*(rhoi/R)**2.
   return factor

#Function for reading and plotting data set.
def readAndPlotFlux(ionData, fluxType, normType):
   #Read flux data from files - see data sets below for indices.
   ionNames        = ionData[names]
   ionFiles        = ionData[files]
   ionStartTs      = ionData[startT]
   ion_qmRats      = ionData[qmRat]
   ion_omts        = ionData[omt]
   ionDensRats     = ionData[densRat]
   ionDataSet      = ionData[dataSet]
   ionDataSetColor = ionData[dataSetColor]
   ionDataFitColor = ionData[dataFitColor]
   avgFlux    = []
   stdDev     = []

   #Read flux data.
   for i,ion in enumerate(ionNames):
      allFluxData = []
      readFlux(ionFiles[i], allFluxData)
      [t, fluxData] = getData(fluxType, allFluxData)
      startIndex = math.ceil(ionStartTs[i]*len(fluxData)) #Get index closest to startT fraction of total time.
      fluxData = np.array(fluxData[startIndex:]) #Need to make a numpy array to perform multiplication by a scalar...
      if (normType == experimentalNorm):
         fluxData = fluxData*getExpFactor()
      elif (normType == impurityNorm):
         fluxData = fluxData * ((1/ion_omts[i])**2) * (1/ionDensRats[i]) #Normalize flux for L_T and n_imp.
      avgFlux.append(getAverageVal(fluxData))
      stdDev.append(getStdDev(fluxData))
   
   #Plot ion points.
   plt.scatter(ion_qmRats, avgFlux, color=ionDataSetColor)
   plt.errorbar(ion_qmRats, avgFlux, yerr=stdDev, linestyle="None", color=ionDataSetColor)

   for i, txt in enumerate(ionNames):
       ax.annotate(txt, (ion_qmRats[i], avgFlux[i] - stdDev[i])) #Offset y for readability.

   #Generate line of best fit.
   ion_qmRatsNew = np.linspace(ion_qmRats[len(ion_qmRats)-1]/1.5, 1, 100) #Note: Extend low end a little past last point.
   coeffs = poly.polyfit(ion_qmRats, avgFlux, 1)
   ffit   = poly.polyval(ion_qmRatsNew, coeffs)
   plt.plot(ion_qmRatsNew, ffit, color=ionDataFitColor, label=ionDataSet)

   textSize = 12
   plt.xlabel("[q/m] / [q/m]$_{H^+}$", fontsize=textSize)
   if (normType == experimentalNorm):
      plt.ylabel("$v_{pinch}$ (m/s)",  fontsize=textSize)
   else:
      plt.ylabel("<$\\Gamma$$_{ES}$>", fontsize=textSize)

#Main routine starts.
normType = ""
if (impNormFlag):
   normType = impurityNorm
elif (expNormFlag):
   normType = experimentalNorm

#Create figure.
fig, ax = plt.subplots()
if (kineticData or bothData):
   readAndPlotFlux(kineticDataSet, fluxType, normType)
if (adiabaticData or bothData):
   readAndPlotFlux(adiabaticDataSet, fluxType, normType)

label  = ""
legend = False
if (kineticData):
   label = "ITG Impurity Pinch - Kinetic Electrons"
   plt.title(label)
elif (adiabaticData):
   label = "ITG Impurity Pinch - Adiabatic Electrons"
   plt.title(label)


plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.tight_layout()
if (kineticData):
   plt.savefig('./GoerlerImpurities/FluxPlotKinEl'   + normType + '.pdf')
elif (adiabaticData):
   plt.savefig('./GoerlerImpurities/FluxPlotAdEl'    + normType + '.pdf')
elif (bothData):
   plt.legend(loc='upper left')
   plt.savefig('./GoerlerImpurities/FluxPlotAdKinEl' + normType + '.pdf')

plt.show()