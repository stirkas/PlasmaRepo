#!/usr/bin/env python

import numpy as np
import numpy.polynomial.polynomial as poly
import math
import matplotlib.pyplot as plt
from plotFlux import readFlux
from plotFlux import getAverageVal
from plotFlux import getData
from plotFlux import getStdDev

dataType  = 0 #Read particle flux.
avgFlux   = []
stdDev    = []
ions      = []
dataFiles = []
qmRat     = []
startT    = []
ltRat     = []
densRat   = []

kineticData = True
if (kineticData):
   ions      = ["H$^+$", "Custom$^{+6}$", "Be$^{+4}$", "Ne$^{+10}$", "Ar$^{+15}$", "Mo$^{+31}$", "W$^{+40}$"]
   dataFiles = ["./GoerlerImpurities/nrgsummary_proton.dat",     "./GoerlerImpurities/nrgsummary_custom.dat",
                "./GoerlerImpurities/nrgsummary_beryllium.dat",  "./GoerlerImpurities/nrgsummary_neon.dat",
                "./GoerlerImpurities/nrgsummary_argon.dat",      "./GoerlerImpurities/nrgsummary_molybdenum.dat",
                "./GoerlerImpurities/nrgsummary_tungsten.dat"]
   qmRat     = [1/1, 6/8, 4/9.0121820, 10/20.17970, 15/39.948, 31/95.95, 40/183.84]
   startT    = [.53, .64, .37, .16, .16, .13, .16] #% location to start at. For ignoring linear phase.
   ltRat     = [6.96, 6.96, 6.96, 6.96, 6.96, 6.96, 6.96] #Ratio of R/L_T = omt to normalize flux different from GENE.
   densRat   = [.001, .001, .001, .001, .001, .001, .001] #Ratio of n_imp/n_e to normalize flux to impurities unlike e- from GENE.
else:
   ions = ["H$^+$", "Ne$^{+10}$", "W$^{+40}$"]
   dataFiles = ["./GoerlerImpurities/nrgsummary_proton_ad.dat", "./GoerlerImpurities/nrgsummary_neon_ad.dat",
                "./GoerlerImpurities/nrgsummary_tungsten_ad.dat"]
   qmRat   = [1/1, 10/20.17970, 40/183.84]
   startT  = [.114, .25, .118]
   ltRat   = [6.96, 6.96, 6.96]
   densRat = [.001, .001, .001]

#Normalization flags. Set at most one true. None gives GENE norm.
impurityNormalization = False  #Convert GENE norm to nicer imp. norm.
experimentalResults   = False  #Convert GENE norm to experimental norm.
normType = ""
if (impurityNormalization):
   normType = "ImpNorm"
elif (experimentalResults):
   normType = "ExpNorm"
experimentalFactor  = []
if (experimentalResults):
   #Begin SEP added stuff.
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

   experimentalFactor = vt/c*(rhoi/R)**2.

   print("rhoi=",rhoi)
   print("vt=",vt)
   print("B=",B)
   print("factor=",experimentalFactor)
   #End of SEP added stuff

for i,ion in enumerate(ions):
   data = []
   readFlux(dataFiles[i], data)
   [t, fluxData] = getData(dataType, data)
   startIndex = math.ceil(startT[i]*len(fluxData)) #Get index closest to startT fraction of total time.
   fluxData = fluxData[startIndex:]
   fluxData = np.array(fluxData) #Need to make a numpy array to perform multiplication by a scalar...
   if (experimentalResults):
      fluxData = fluxData*experimentalFactor
   elif (impurityNormalization):
      fluxData = fluxData * ((1/ltRat[i])**2) * (1/densRat[i]) #Normalize flux for L_T and n_imp.
   avgFlux.append(getAverageVal(fluxData))
   stdDev.append(getStdDev(fluxData))

#Plot ion points and best fit line.
fig, ax = plt.subplots()
plt.scatter(qmRat, avgFlux)
plt.errorbar(qmRat, avgFlux, yerr=stdDev, linestyle="None")

for i, txt in enumerate(ions):
    ax.annotate(txt, (qmRat[i], avgFlux[i] - stdDev[i])) #Offset y for readability.

#Generate line of best fit.
qmRatNew = np.linspace(qmRat[len(qmRat)-1]/1.5, 1, 100) #Note: Extend low end a little past last point.
coeffs = poly.polyfit(qmRat, avgFlux, 1)
ffit   = poly.polyval(qmRatNew, coeffs)
label = ""
if (kineticData):
   label = "Impurity Pinch - ITG w/ Kinetic Electrons"
else:
   label = "Impurity Pinch - ITG w/ Adiabatic Electrons"
plt.title(label)
plt.plot(qmRatNew, ffit, color='r')

plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
textSize = 12
plt.xlabel("[q/m] / [q/m]$_{H^+}$", fontsize=textSize)
if (experimentalResults):
   plt.ylabel("$v_{pinch}$ (m/s)",  fontsize=textSize)
else:
   plt.ylabel("<$\\Gamma$$_{ES}$>", fontsize=textSize)
plt.tight_layout()
if (kineticData):
   plt.savefig('./GoerlerImpurities/FluxPlotKinEl' + normType + '.pdf')
else:
   plt.savefig('./GoerlerImpurities/FluxPlotAdEl'  + normType + '.pdf')
plt.show()