#!/usr/bin/env python

#Numerical calculation of flux from Estrada (05) - eq.'s 14 & 15.
#Not sure what to use for phi...going to take phi, s-hat, etc. from GENE ETG run to match flux data.

import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
import shutil
import sys
from scipy.integrate import dblquad
import time

#Do this to import from geneScripts
sys.path.append('../geneScripts')
from plotGrowthRates import readGrowthRates

#Define a function to integrate real and imag parts separately.
def complex_quadrature(func, xa, xb, ya, yb, **kwargs):
   def real_func(y,x):
      return scipy.real(func(y,x))
   def imag_func(y,x):
      return scipy.imag(func(y,x))
   real_integral = dblquad(real_func, xa, xb, ya, yb, **kwargs)
   imag_integral = dblquad(imag_func, xa, xb, ya, yb, **kwargs)
   return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

adiabaticPureFile     = '../geneScripts/GoerlerImpLin/scanfiles0000/scan.log'
adiabaticProtImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0001/scan.log'
adiabaticNeonImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0002/scan.log'
adiabaticTungImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0003/scan.log'
adiabaticPure     = readGrowthRates(adiabaticPureFile)
adiabaticProtImp  = readGrowthRates(adiabaticProtImpFile)
adiabaticNeonImp  = readGrowthRates(adiabaticNeonImpFile)
adiabaticTungImp  = readGrowthRates(adiabaticTungImpFile)
linearGrowthData = [adiabaticPure, adiabaticProtImp, adiabaticNeonImp, adiabaticTungImp]

#Geometric data.
R     = 1    #Major radius of tokamak.
Lref  = R    #To match GENE with omn/omt.
s     = .837 #Shear parameter - s-hat
#Shared data for allocating data arrays.
L_T      = Lref/6.96 #L_T for ions and electrons - Lref/omt.
zI       = [1, 10, 40]
mI       = [1, 20.179700, 183.84]
nI       = .001 #Impurity density / elec density.
m_i      = 1
m_e      = 5.4551*10**-4 #m_e / m_i technically.
z_e      = -1
z_i      = -z_e

#Load data
#Species indices: i = main ion, e = electron, I = impurity
#Data indices: 0 = z, 1 = m (amu), 2 = n (n_e), 3 = L_T, 4 = L_n
#f = z*n/n_e, mu = sqrt(m_i/m_I)
#Data for each main species
ionDensities   = [.999, .99, .96]
ionDensGrads   = [Lref/2.2222222, Lref/2.424242, Lref/2.3125] #L_n for imps - Lref/omn.
mainIonData    = [z_i,  m_i,  ionDensities,   L_T, ionDensGrads, ionDensities,     1] #Protons
electronData   = [z_e,  m_e,  1,              L_T, Lref/2.22,    1, np.sqrt(m_i/m_e)] #Note: L_T,e always equal to L_T,i according to this theory.

#Impurity data.
protonData     = [zI[0], mI[0], nI, np.inf, np.inf]
neonData       = [zI[1], mI[1], nI, np.inf, np.inf]
tungstenData   = [zI[2], mI[2], nI, np.inf, np.inf]
#Store all data together, put ions and e at the end.
allData = [protonData, neonData, tungstenData, mainIonData, electronData]
speciesNames = ["proton", "neon", "tungsten", "main ion", "electrons"]
numImpurities = allData - 2
flux   = np.zeros((numImpurities, len(kthetaRhoi)))
flux_i = np.zeros((numImpurities, len(kthetaRhoi)))

phi   = 50/474 #Stolen from GENE ITG imp. cases. Multiplied to go to general normalization not GENE's rho* version.

#Average over theta. Dependence only in curvature/gradient drift velocity term. 1/2pi already in moved to other terms.
thetaFunc     = lambda x: np.cos(x) + s*x*np.sin(x)
thetaIntegral = quad(thetaFunc, 0, 2*np.pi)

cgsUnits = True
if (cgsUnits):
   m_i = m_i * 1.6726219*10**-24 #grams
   mI  = mI  * m_i
   me  = m_e * m_i
   z_e = z_e * 4.8*10**-10 #cgs, stat-couloumbs
   z_i = -z_e
   zI = zI * z_i
   T_e = 1.6*10**-9 #ergs - 1keV
   T_i = T_e
   T_I = .1*T_e
   B = 10000 #Gauss - 1 Tesla

for i, species in enumerate(speciesNames):
   #Gather response function data for each impurity and the main ion.
   zI   = allData[i][0]
   L_TI = allData[i][3]
   LnI  = allData[i][4]
   fI   = allData[i][5]
   muI  = allData[i][6]

   if (i == numImpurities): #Kick out early once we're at ions.
      break

   z_i  = allData[numImpurities][0]
   Ln_i = allData[numImpurities][4][i]
   f_i  = allData[numImpurities][5][i]
   mu_i = allData[numImpurities][6]
   
   print('imp: ' + speciesNames[i])

#   # y = vPerp, x = vPar - will evaluate y integral first.
#   for j, kthRho in enumerate(kthetaRhoi):
#      omega        = omegaReal[j] + 1j*gamma[j]
#      oStarI       = lambda y, x: kthRho * (1/LnI  + (y**2 + x**2 - 3)/(L_TI*2)) #omega_star,I - Eq. (10)
#      oStar_i      = lambda y, x: kthRho * (1/Ln_i + (y**2 + x**2 - 3)/(L_T*2))
#      oD_I         = lambda y, x: (2*kthRho/(zI*R))  * ((y**2)/4 + (x**2)/2) * thetaIntegral #omega_d,I - Eq. (11)
#      oDi          = lambda y, x: (2*kthRho/(z_i*R)) * ((y**2)/4 + (x**2)/2) * thetaIntegral
#      omegaTerm    = lambda y, x: (zI*omega  + oStarI(y,x))/(omega  + oD_I(y,x))
#      omegaTerm_i  = lambda y, x: (z_i*omega + oStar_i(y,x))/(omega + oDi(y,x))
#      xTerm        = lambda x: scipy.exp(-1*(x**2)/2)   #vPar  int.
#      yTerm        = lambda y: y*scipy.exp(-1*(y**2)/2) #vPerp int.
#      besselTerm   = lambda y: scipy.special.jv(0, kthRho*y/(zI*muI))**2
#      besselTerm_i = lambda y: scipy.special.jv(0, kthRho*y/(z_i*mu_i))**2
#
#      func     = lambda y, x: xTerm(x)*yTerm(y)*besselTerm(y)*omegaTerm(y,x)
#      result   = complex_quadrature(func, -scipy.inf, scipy.inf, 0, scipy.inf)
#      func_i   = lambda y, x: xTerm(x)*yTerm(y)*besselTerm_i(y)*omegaTerm_i(y,x)
#      result_i = complex_quadrature(func_i, -scipy.inf, scipy.inf, 0, scipy.inf)
#      
#      #Response function for impurities.
#      R_I[i][j]  = -zI*fI + (fI/np.sqrt(2*np.pi))*result[0]
#      flux[i][j] = np.real(R_I[i][j] * kthRho * phi**2 * 1j)
#      R_i[i][j]  = -z_i*f_i + (f_i/np.sqrt(2*np.pi))*result_i[0]
#      flux_i[i][j] = np.real(R_i[i][j] * kthRho * phi**2 * 1j)
#      print('i: '     + str(j))
#      print(flux[i][j])
#      print(flux_i[i][j])


dirName = 'GoerlerImpurities'
if not os.path.exists(dirName):
    os.makedirs(dirName)
else:
   shutil.rmtree(dirName)
   os.makedirs(dirName)

np.savez_compressed(dirName + '/flux.npz',       flux)
np.savez_compressed(dirName + '/species.npz',    speciesNames[:-2])
np.savez_compressed(dirName + '/growthData.npz', linearGrowthData)