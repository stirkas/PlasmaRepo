#!/usr/bin/env python

import numpy as np
import numpy.polynomial.polynomial as poly
import math
import matplotlib.pyplot as plt
from plotFlux import readFlux
from plotFlux import getFluxAverage
from plotFlux import getFluxData
from plotFlux import getFluxStdDev

#Set data flags for script.
#User input flags - just change these input flags if you want to ignore every choice until runtime.
allowGeneralInput = False #Accept input for all other flags.
allowParamInput   = False #Accept tokamak param. input.
#Data type flag. Set only one.
adiabaticData = True
kineticData   = False
#Normalization flags. Set at most one true. None gives GENE norm.
impNormFlag = False #Convert GENE norm to nicer impurity norm using L_T and n_imp.
expNormFlag = False #Convert GENE norm to experimental norm.
#Flux data types
particleFlux = 0
heatFlux     = 1
momentumFlux = 2
fluxType  = particleFlux

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
invertX      = 9
invertY      = 10
title        = 11
species = ['Bulk Ions', 'Electrons', 'Impurities']
#Normalization types
experimentalNorm = "ExpNorm"
impurityNorm     = "ImpNorm"
geneNorm         = "GB" #GyroBohm norm.

#Data sets.
kineticIonData =   [["H$^+$", "Custom$^{+6}$", "Be$^{+4}$", "Ne$^{+10}$", "Ar$^{+15}$", "Mo$^{+31}$", "W$^{+40}$"],
                    ["./GoerlerImpurities/nrgsummary_protonIon.dat",     "./GoerlerImpurities/nrgsummary_customIon.dat",
                     "./GoerlerImpurities/nrgsummary_berylliumIon.dat",  "./GoerlerImpurities/nrgsummary_neonIon.dat",
                     "./GoerlerImpurities/nrgsummary_argonIon.dat",      "./GoerlerImpurities/nrgsummary_molybdenumIon.dat",
                     "./GoerlerImpurities/nrgsummary_tungstenIon.dat"],
                    [1/1, 6/8, 4/9.0121820, 10/20.17970, 15/39.948, 31/95.95, 40/183.84],
                    [.001, .001, .001, .001, .001, .001, .001],
                    [.53, .64, .37, .16, .16, .13, .16],
                    [6.96, 6.96, 6.96, 6.96, 6.96, 6.96, 6.96],
                    "kinetic",
                    "steelblue",
                    "red",
                    True,
                    False,
                    species[0]]
kineticElData =    [["H$^+$", "Custom$^{+6}$", "Be$^{+4}$", "Ne$^{+10}$", "Ar$^{+15}$", "Mo$^{+31}$", "W$^{+40}$"],
                    ["./GoerlerImpurities/nrgsummary_protonEl.dat",     "./GoerlerImpurities/nrgsummary_customEl.dat",
                     "./GoerlerImpurities/nrgsummary_berylliumEl.dat",  "./GoerlerImpurities/nrgsummary_neonEl.dat",
                     "./GoerlerImpurities/nrgsummary_argonEl.dat",      "./GoerlerImpurities/nrgsummary_molybdenumEl.dat",
                     "./GoerlerImpurities/nrgsummary_tungstenEl.dat"],
                    [1/1, 6/8, 4/9.0121820, 10/20.17970, 15/39.948, 31/95.95, 40/183.84],
                    [.001, .001, .001, .001, .001, .001, .001],
                    [.53, .64, .37, .16, .16, .13, .16],
                    [6.96, 6.96, 6.96, 6.96, 6.96, 6.96, 6.96],
                    "kinetic",
                    "steelblue",
                    "red",
                    True,
                    False,
                    species[1]]
kineticImpData =   [["H$^+$", "Custom$^{+6}$", "Be$^{+4}$", "Ne$^{+10}$", "Ar$^{+15}$", "Mo$^{+31}$", "W$^{+40}$"],
                    ["./GoerlerImpurities/nrgsummary_protonImp.dat",     "./GoerlerImpurities/nrgsummary_customImp.dat",
                     "./GoerlerImpurities/nrgsummary_berylliumImp.dat",  "./GoerlerImpurities/nrgsummary_neonImp.dat",
                     "./GoerlerImpurities/nrgsummary_argonImp.dat",      "./GoerlerImpurities/nrgsummary_molybdenumImp.dat",
                     "./GoerlerImpurities/nrgsummary_tungstenImp.dat"],
                    [1/1, 6/8, 4/9.0121820, 10/20.17970, 15/39.948, 31/95.95, 40/183.84],
                    [.001, .001, .001, .001, .001, .001, .001],
                    [.53, .64, .37, .16, .16, .13, .16],
                    [6.96, 6.96, 6.96, 6.96, 6.96, 6.96, 6.96],
                    "kinetic",
                    "steelblue",
                    "red",
                    True,
                    True,
                    species[2]]
kineticDataSet   = [kineticIonData, kineticElData, kineticImpData]
adiabaticIonData = [["H$^+$", "Ne$^{+10}$", "W$^{+40}$"],
                    ["./GoerlerImpurities/nrgsummary_protonIon_ad.dat", "./GoerlerImpurities/nrgsummary_neonIon_ad.dat",
                     "./GoerlerImpurities/nrgsummary_tungstenIon_ad.dat"],
                    [1/1, 10/20.17970, 40/183.84],
                    [.001, .001, .001],
                    [.114, .143, .125],
                    [6.96, 6.96, 6.96],
                    "adiabatic",
                    "red",
                    "blue",
                    True,
                    False,
                    species[0]]
adiabaticImpData = [["H$^+$", "Ne$^{+10}$", "W$^{+40}$"],
                    ["./GoerlerImpurities/nrgsummary_protonImp_ad.dat", "./GoerlerImpurities/nrgsummary_neonImp_ad.dat",
                     "./GoerlerImpurities/nrgsummary_tungstenImp_ad.dat"],
                    [1/1, 10/20.17970, 40/183.84],
                    [.001, .001, .001],
                    [.114, .143, .125],
                    [6.96, 6.96, 6.96],
                    "adiabatic",
                    "red",
                    "blue",
                    True,
                    True,
                    species[2]]
adiabaticDataSet = [adiabaticIonData, adiabaticImpData]

#Accept user input before running main routine.
# Parameters are all obtained from Callen Nuc. Fusion 2010 by Scott
# This is a DIII-D H-mode before the onset of an ELM
R0  = 1.7    # Magnetic axis, meters
R   = 2.2    # location of cold impurities, near separatrix
B0  = 2.     # Teslas, B at R0
mu  = 2.     # m/mp for bulk ions - I use mu=2 because tokamaks typically run with D.
Z   = 1.     # This is Z for the bulk ions
Ti  = 1000.  # bulk ion temp in eV at cold impurity location
c   = 0.001  # c=nI/n0, concentration of impurity
if (allowGeneralInput):
   #Reset all flags so defaults above don't screw things up.
   adiabaticData = False
   kineticData   = False
   impNormFlag   = False
   expNormFlag   = False
   #Read flags from input.
   adiabaticData = input("Plot data for adiabatic electrons? (y/n): ") == "y"
   kineticData   = input("Plot data for kinetic electrons? (y/n): ")   == "y"
   dataType      = int(input("Choose a flux type. (0 = particle flux (only one working), 1 = heat flux, 2 = momentum flux): "))
   normalization = int(input("Choose a normalization. (0 = GENE, 1 = impurity (L_ref = L_T, n_ref = n_I), 2 = experimental): "))
   #Reset normalization flags and set from input.
   if (normalization == 1):
      impNormFlag = True
   elif (normalization == 2):
      expNormFlag = True
   

if (allowParamInput):
   paramInput = input("Would you like to input necessary tokamak parameters? (y/n): ")
   if (paramInput == "y"):
      R0 = float(input("Major radius R0 (m): "))
      R  = float(input("Edge radius R (m): "))
      B0 = float(input("Background B field (T): "))
      mu = float(input("Bulk ion to proton mass ratio: "))
      Z  = float(input("Bulk ion charge factor Z: "))
      Ti = float(input("Bulk ion temp (eV): "))
      c  = float(input("Impurity concentration, n_I/n_e: "))

#Function for converting flux to SI units with DIII-D params - by Scott.
def getExpFactor():
   # Flux = Flux_GENE x n0 vt (rhoi/R)^2
   # R=L_ref, n0=n_ref, vt=vt_ref, rhoi_ref
   # vt=979.0/sqrt(mu) x sqrt(Ti)
   # rhoi=1.02 x sqrt(mu Ti)/(Z B)
   # factor= Flux/(nI*Flux_GENE) = 1/c vt (rhoi/R)^2

   # a little tricky... We want v_pinch.  Flux = v_pinch x nI
   # Flux_GENE = v_pinch x n0.  v_pinch = Flux_GENE/(nI/n0) x vt (rhoi/R)^2

   B = B0*R0/R*10000. # Plasma formulary formulas use Gauss
   rhoi = 1.02*np.sqrt(mu*Ti)/(Z*B)
   vti=9790.*np.sqrt(Ti/mu)

   factor = vti/c*(rhoi/R)**2.
   return factor

#Function for reading and plotting data set.
def readAndPlotFlux(ionData, axis):
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
   invertDataX     = ionData[invertX]
   invertDataY     = ionData[invertY]
   dataTitle       = ionData[title]
   avgFlux         = []
   stdDev          = []

   #Read flux data.
   for i,ion in enumerate(ionNames):
      allFluxData = []
      readFlux(ionFiles[i], allFluxData)
      [t, fluxData] = getFluxData(fluxType, allFluxData)
      startIndex = math.ceil(ionStartTs[i]*len(fluxData)) #Get index closest to startT fraction of total time.
      fluxData = np.array(fluxData[startIndex:]) #Need to make a numpy array to perform multiplication by a scalar...
      if (expNormFlag):
         fluxData = fluxData*getExpFactor()
      elif (impNormFlag):
         fluxData = fluxData * ((1/ion_omts[i])**2) * (1/ionDensRats[i]) #Normalize flux for L_T and n_imp.
      avgFlux.append(getFluxAverage(fluxData))
      stdDev.append(getFluxStdDev(fluxData))
   
   #Plot ion points.
   axis.scatter(ion_qmRats, avgFlux, color=ionDataSetColor)
   axis.errorbar(ion_qmRats, avgFlux, yerr=stdDev, linestyle="None", color=ionDataSetColor)

   for i, txt in enumerate(ionNames):
       axis.annotate(txt, (ion_qmRats[i], avgFlux[i] - stdDev[i])) #Offset y for readability.

   #Generate line of best fit.
   ion_qmRatsNew = np.linspace(ion_qmRats[-1]/1.5, 1, 100) #Note: Extend low end a little past last point.
   coeffs = poly.polyfit(ion_qmRats, avgFlux, 1)
   ffit   = poly.polyval(ion_qmRatsNew, coeffs)
   axis.plot(ion_qmRatsNew, ffit, color=ionDataFitColor, label=ionDataSet)

   textSize = 12
   largeTextSize = 16
   axis.set_xlabel("[q/m] / [q/m]$_{H^+}$", fontsize=textSize)
   if (expNormFlag):
      axis.set_ylabel("$v_{pinch}$ (m/s)",  fontsize=textSize)
   elif (impNormFlag):
      axis.set_ylabel("$\\frac{<\\Gamma>}{n_Iv_{Ti}(\\rho_i/L_{Ti})^2}$", fontsize=largeTextSize)
   else:
      axis.set_ylabel("<$\\Gamma$$_{GB}$>", fontsize=textSize)
   
   if (invertDataX):
      axis.invert_xaxis()
   if (invertDataY):
      axis.invert_yaxis()

   axis.set_title(dataTitle)
   axis.grid()

#Main routine starts.

#Create figure and plot data.
if (kineticData):
   fig, axs = plt.subplots(1, len(kineticDataSet))
   fig.set_figheight(6)
   fig.set_figwidth(12)
   for i in range(len(kineticDataSet)):
      readAndPlotFlux(kineticDataSet[i], axs[i])
   axs[1].set_xlim(axs[0].get_xlim())
   axs[1].set_ylim(axs[0].get_ylim())
elif (adiabaticData):
   fig, axs = plt.subplots(1, len(adiabaticDataSet))
   for i in range(len(adiabaticDataSet)):
      readAndPlotFlux(adiabaticDataSet[i], axs[i])

label  = ""
if (adiabaticData):
   fig.suptitle("Electrostatic ITG Impurity Pinch - Adiabatic Electrons\n[$n_I=.001n_e, T_I=.1T_i=.1T_e$]")
elif (kineticData):
   fig.suptitle("Electrostatic ITG Impurity Pinch - Kinetic Electrons\n[$n_I=.001n_e, T_I=.1T_i=.1T_e$]")

fig.tight_layout(rect=[0, 0.03, 1, 0.93])

#Get norm type in string form for file saving.
normType = ""
if (impNormFlag):
   normType = impurityNorm
elif (expNormFlag):
   normType = experimentalNorm
else:
   normType = geneNorm

if (adiabaticData):
   plt.savefig('./GoerlerImpurities/FluxPlotAdEl'    + normType + '.pdf')
elif (kineticData):
   plt.savefig('./GoerlerImpurities/FluxPlotKinEl'   + normType + '.pdf')

plt.show()