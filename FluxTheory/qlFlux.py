#!/usr/bin/env python

#Numerical calculation of flux from Estrada (05) - eq.'s 14 & 15.
#Not sure what to use for phi...going to take phi, s-hat, etc. from GENE ETG run to match flux data.

import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
import shutil
import sys
from scipy.integrate import quad
from scipy.integrate import dblquad
import time

#TODO: Add ballooning/ambipolairty to QL theory.

#Do this to import from geneScripts
sys.path.append('../geneScripts')
from plotGrowthRates import readGrowthRates
from readAmpSpectra  import readAmpSpectra

#Define a function to integrate real and imag parts separately.
def complex_quadrature(func, xa, xb, ya, yb, **kwargs):
   def real_func(y,x):
      return scipy.real(func(y,x))
   def imag_func(y,x):
      return scipy.imag(func(y,x))
   real_integral = dblquad(real_func, xa, xb, ya, yb, **kwargs)
   imag_integral = dblquad(imag_func, xa, xb, ya, yb, **kwargs)
   return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

adiabaticProtImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0001/scan.log'
adiabaticNeonImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0002/scan.log'
adiabaticTungImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0003/scan.log'
adiabaticProtImp  = readGrowthRates(adiabaticProtImpFile)
adiabaticNeonImp  = readGrowthRates(adiabaticNeonImpFile)
adiabaticTungImp  = readGrowthRates(adiabaticTungImpFile)

ampSpectraFiles = ['../geneScripts/GoerlerImpurities/AmpSpectraProtonImpAdEl.dat',
                   '../geneScripts/GoerlerImpurities/AmpSpectraNeonImpAdEl.dat',
                   '../geneScripts/GoerlerImpurities/AmpSpectraTungstenImpAdEl.dat']

#Choose between ad or kin electron data.
adiabatic = True
numModes  = 0
if (adiabatic):
   linearGrowthData = [adiabaticProtImp, adiabaticNeonImp, adiabaticTungImp]
else:
   linearGrowthData = 0 #Read kinetic data when ready.
numModes = numModes - 1 #Remove last mode because phi mode data starts from krho = 0

#Geometric data.
s        = .837 #Shear parameter - s-hat
#Shared data for allocating data arrays.
n_e   = 1 #Normalize to n_e = 1.
nI    = .001*n_e #Impurity density / elec density really.
m_i   = 1
m_e   = 5.4462*10**-4 #m_e / m_i really. Same # used in GENE params.
z_e   = -1
z_i   = -z_e
Ti    = 1    #Ion temperature = reference temp.
Te    = Ti 
T_I   = .1*Ti
vTi   = np.sqrt(Ti/m_i)
zIall = np.array([1, 10, 40], dtype = float)          #Norm to main ion charge.
mIall = np.array([1, 20.179700, 183.84], dtype=float) #Norm to main ion mass.

#Load data
#Species indices: i = main ion, e = electron, I = impurity
#Data indices: 0 = z, 1 = m (amu), 2 = n (n_e)
#Data for each main species
ionDensities   = np.array([.999, .99, .96])*n_e
omn_i          = np.array([2.2222222, 2.2424242, 2.3125])
omt_i          = 6.96
omnI           = 0
omtI           = 0
mainIonData    = [z_i, m_i, ionDensities] #Protons for main ion.
electronData   = [z_e, m_e,            1]

#Impurity data.
protonData     = [zIall[0], mIall[0], nI]
neonData       = [zIall[1], mIall[1], nI]
tungstenData   = [zIall[2], mIall[2], nI]
#Store all data together, put ions and e at the end.
allData = [protonData, neonData, tungstenData, mainIonData, electronData]
speciesNames = ["proton", "neon", "tungsten", "main ion", "electrons"]
numImpurities = len(allData) - 2 #2 for main ion species and electrons.

flux   = np.zeros(numImpurities, dtype=complex)
flux_i = np.zeros(numImpurities, dtype=complex)

eta = 0 #Ballooning angle. Phi peaks at 0 so use 0 for now.
#Peaked phi at 0 eta for each mode for each impurity.
#phi = [[.8, .45,  .125, .225, .25, .2,   .125, .04, .05, .04, .02, .04, .02,  .04, .06],
#       [.4, .225, .45,  .15,  .15, .125, .07,  .2,  .04, .04, .04, .04, .05,  .05, .015],
#       [.1, .285, .12,  .25,  .06, .2,   .08,  .07, .04, .05, .04, .02, .015, .02, .015]]

#If not using ballooning spectra, integrate theta dependence.
thetaFunc = lambda x: (1/(2*np.pi))*(np.cos(x) - s*x*np.sin(x))
thetaInt  = quad(thetaFunc, 0, 2*np.pi)[0]

ampSpectra = []

for i, species in enumerate(speciesNames):
   if (i == numImpurities): #Kick out early once we're passed all impurities.
      break

   print('------')
   print(species)

   kthetaRhoi = linearGrowthData[i][0]
   omegaReal  = linearGrowthData[i][1]
   gamma      = linearGrowthData[i][2]

   #Gather response function data for each impurity and the main ion.
   zI   = allData[i][0]
   mI   = allData[i][1]
   nI   = allData[i][2]
   vT_I = np.sqrt(T_I/mI)
   cycI = zI/(mI) #Cyc freq. #Dont need B or c because cyc freq always shows up in ratios of main ion cyc freq and consts will cancel.
   rhoI = vT_I/cycI

   cyc_i = (z_i)/(m_i) #Cyc freq. #Dont need B or c because cyc freq always shows up in ratios of main ion cyc freq and consts will cancel.
   rho_i = vTi/cyc_i
   n_i   = ionDensities[i]

   #Phi for ballooning angle spectra.
   #phiImp = phi[i] #Drop krho = 0 term.
   
   #Load phi from GENE here for the specific run. Then split by mode.
   ampSpectra = readAmpSpectra(ampSpectraFiles[i])
   phiAmpImp  = np.array(ampSpectra[1][1])
   ampSpectra.append(phiAmpImp[1:])

   #Loop over phi since it is missing the highest kthetaRho term so will not break on last term.
   #y = vPerp, x = vPar - will evaluate y integral first.
   for j, phiMode in enumerate(phiAmpImp):
      omega        = (omegaReal[j] + 1j*gamma[j])
      oStarT_I     = lambda y, x: -kthetaRhoi[j] * (vT_I/vTi) * (rhoI/rho_i)  * (omnI     + (1/2)*((y**2)*((vTi/vT_I)**2) + (x**2)*((vTi/vT_I)**2) - 3)*omtI)
      oStarT_i     = lambda y, x: -kthetaRhoi[j] * (vTi/vTi)  * (rho_i/rho_i) * (omn_i[i] + (1/2)*((y**2)*((vTi/vTi)**2)  + (x**2)*((vTi/vTi)**2)  - 3)*omt_i)
      obarD_I      = lambda y, x:  kthetaRhoi[j] * (vT_I/vTi) * (rhoI/rho_i)  * ((1/2)*(y**2)  + (x**2)) * ((vTi/vT_I)**2) * thetaInt #(np.cos(eta) + s*eta*np.sin(eta)) For ballooning modes.
      obarDi       = lambda y, x:  kthetaRhoi[j] * (vTi/vTi)  * (rho_i/rho_i) * ((1/2)*(y**2)  + (x**2)) * ((vTi/vTi)**2)  * thetaInt #(np.cos(eta) + s*eta*np.sin(eta)) For ballooning modes.
      omegaTermI   = lambda y, x: (omega - oStarT_I(y,x))/(omega - obarD_I(y,x)) - (np.conj(omega) - oStarT_I(y,x))/(np.conj(omega) - obarD_I(y,x))
      omegaTerm_i  = lambda y, x: (omega - oStarT_i(y,x))/(omega - obarDi(y,x))  - (np.conj(omega) - oStarT_i(y,x))/(np.conj(omega) - obarDi(y,x))
      xTerm_i      = lambda x: x**2*scipy.exp((-1/2)*(x**2)*((vTi/vTi)**2))
      yTerm_i      = lambda y: y**3*scipy.exp((-1/2)*(y**2)*((vTi/vTi)**2))
      xTermI       = lambda x: x**2*scipy.exp((-1/2)*(x**2)*((vTi/vT_I)**2))
      yTermI       = lambda y: y**3*scipy.exp((-1/2)*(y**2)*((vTi/vT_I)**2))
      besselTermI  = lambda y: scipy.special.jv(0, kthetaRhoi[j] * y *  (cyc_i/cycI) * np.sqrt(1 + (s*eta)**2))**2
      besselTerm_i = lambda y: scipy.special.jv(0, kthetaRhoi[j] * y * (cyc_i/cyc_i) * np.sqrt(1 + (s*eta)**2))**2
      constantTermI  = .5j * kthetaRhoi[j] * zI  * (Ti/T_I) * (mI/m_i)  * nI  * phiMode**2 * (1/(2*np.pi))**(1/2) * (vTi/vT_I)**3 #phiImp[j] for ballooning modes.
      constantTerm_i = .5j * kthetaRhoi[j] * z_i * (Ti/Ti)  * (m_i/m_i) * n_i * phiMode**2 * (1/(2*np.pi))**(1/2) * (vTi/vTi)**3  #phiImp[j] for ballooning modes.
      func     = lambda y, x: constantTermI*xTermI(x)*yTermI(y)*besselTermI(y)*omegaTermI(y,x)
      result   = complex_quadrature(func, -scipy.inf, scipy.inf, 0, scipy.inf)[0]
      func_i   = lambda y, x: constantTerm_i*xTerm_i(x)*yTerm_i(y)*besselTerm_i(y)*omegaTerm_i(y,x)
      result_i = complex_quadrature(func_i, -scipy.inf, scipy.inf, 0, scipy.inf)[0]
      flux[i]   = flux[i]   + result
      flux_i[i] = flux_i[i] + result_i
      print(kthetaRhoi[j])

#Normalize ballooning peak in terms of amp spectra peak.
#ampSpectra = [12,9,15]
for i, species in enumerate(speciesNames):
   if (i==numImpurities):
      break
   print('-------')
   print(species)
   print(flux[i])
   print(flux_i[i])
   #print(flux[i]*((15/phi[2][j])**2))
   #print(np.real(flux_i[i])*((15/phi[2][j])**2))

#Save data off.
#dirName = 'GoerlerImpurities'
#if not os.path.exists(dirName):
#    os.makedirs(dirName)
#else:nc_i, -scipy.inf
#   shutil.rmtree(dirName)
#   os.makedirs(dirName)
#
#np.savez_compressed(dirName + '/flux.npz',       flux)
#np.savez_compressed(dirName + '/flux_i.npz',     flux_i)
#np.savez_compressed(dirName + '/impSpecies.npz', speciesNames[:-2])
#np.savez_compressed(dirName + '/growthData.npz', linearGrowthData)