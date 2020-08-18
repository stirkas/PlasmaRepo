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

adiabaticPureFile     = '../geneScripts/GoerlerImpLin/scanfiles0000/scan.log'
adiabaticProtImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0001/scan.log'
adiabaticNeonImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0002/scan.log'
adiabaticTungImpFile  = '../geneScripts/GoerlerImpLin/scanfiles0003/scan.log'
adiabaticProtImp  = readGrowthRates(adiabaticProtImpFile)
adiabaticNeonImp  = readGrowthRates(adiabaticNeonImpFile)
adiabaticTungImp  = readGrowthRates(adiabaticTungImpFile)
adiabaticPure     = readGrowthRates(adiabaticPureFile)

ampSpectraFiles = ['../geneScripts/GoerlerImpurities/ampSpectraAdProtons.dat',
                   '../geneScripts/GoerlerImpurities/ampSpectraAdNeon.dat',
                   '../geneScripts/GoerlerImpurities/ampSpectraAdTungsten.dat',]

#Choose between ad or kin electron data.
adiabatic = True
numModes  = 0
if (adiabatic):
   numModes = 16
   linearGrowthData = [adiabaticProtImp, adiabaticNeonImp, adiabaticTungImp, adiabaticPure]
else:
   numModes = 23
   linearGrowthData = 0
numModes = numModes - 1 #Remove last mode because phi mode data starts from krho = 0

#Geometric data.
R     = 1    #Major radius of tokamak. Data files say 1.65 though...
s     = .837 #Shear parameter - s-hat
#Shared data for allocating data arrays.
zIall    = np.array([1, 10, 40], dtype = float)
mIall    = np.array([1, 20.179700, 183.84], dtype=float)
n_e      = 10**19   # #/m^3
nI       = .001*n_e #Impurity density / elec density.
m_i      = 1
m_e      = 5.4551*10**-4 #m_e / m_i technically.
z_e      = -1
z_i      = -z_e
c        = 3*10**8 #Speed of light
B        = 1 #Tesla
Te       = 1.6*10**-16 #Joules - 1keV

#Allow for non-normalized units.
cgsUnits = True
if (cgsUnits):
   m_i   = m_i * 1.6726219*10**-24 #grams
   mIall = mIall  * m_i
   me    = m_e * m_i
   z_e   = z_e * 4.8*10**-10 #cgs, stat-couloumbs
   z_i   = -z_e
   zIall = zIall * z_i
   Te    = Te*10**7 #Joules to ergs
   B     = B*10000 #Gauss - 1 Tesla
   c     = c*100
   R     = R*100
   n_e   = n_e*10^-6

Ti   = Te
T_I  = .1*Te
vTi  = np.sqrt(Ti/m_i)
Lref = R         #To match GENE with omn/omt.
L_T  = Lref/6.96 #L_T for ions and electrons - Lref/omt.

#Load data
#Species indices: i = main ion, e = electron, I = impurity
#Data indices: 0 = z, 1 = m (amu), 2 = n (n_e), 3 = L_T, 4 = L_n
#f = z*n/n_e, mu = sqrt(m_i/m_I)
#Data for each main species
ionDensities   = np.array([.999, .99, .96])*n_e
ionLn          = [Lref/2.2222222, Lref/2.424242, Lref/2.3125, Lref/2.22] #L_n for imps - Lref/omn. No imps last.
mainIonData    = [z_i, m_i, ionDensities, L_T,     ionLn] #Protons
electronData   = [z_e, m_e, 1,            L_T, Lref/2.22]

#Impurity data.
protonData     = [zIall[0], mIall[0], nI, np.inf, np.inf]
neonData       = [zIall[1], mIall[1], nI, np.inf, np.inf]
tungstenData   = [zIall[2], mIall[2], nI, np.inf, np.inf]
#Store all data together, put ions and e at the end.
allData = [protonData, neonData, tungstenData, mainIonData, electronData]
speciesNames = ["proton", "neon", "tungsten", "main ion", "electrons"]
numImpurities = len(allData) - 2 #2 for main ion species and electrons.

kthetaRhoi = linearGrowthData[3][0] #Grab lowest index from any data set. kth*rho same for all.
#Technically tung ad longer than rest. But gamma and omega go to 0 for higher ky*rho than the 16th.
flux   = np.zeros((numImpurities, len(kthetaRhoi)))
flux_i = np.zeros((numImpurities, len(kthetaRhoi)))

#Average over theta. Dependence only in curvature/gradient drift velocity term. 1/2pi already in moved to other terms.
thetaFunc     = lambda x: np.cos(x) + s*x*np.sin(x)
thetaIntegral = quad(thetaFunc, 0, 2*np.pi)[0] #0 takes only the real part. Since this integral should just be real.

for i, species in enumerate(speciesNames):
   if (i == numImpurities): #Kick out early once we're passed all impurities.
      break

   kthetaRhoi = linearGrowthData[i][0]
   omegaReal  = linearGrowthData[i][1]
   gamma      = linearGrowthData[i][2]

   #Gather response function data for each impurity and the main ion.
   zI   = allData[i][0]
   mI   = allData[i][1]
   nI   = allData[i][2]
   L_TI = allData[i][3]
   LnI  = allData[i][4]
   vT_I = np.sqrt(T_I/mI)
   cycI = zI*B/(c*mI) #Cyc freq.
   rhoI = vT_I/cycI

   cyc_i = (z_i*B)/(c*m_i) #Cyc freq.
   rho_i = vTi/cyc_i
   Ln_i  = allData[numImpurities][4][i]
   n_i   = ionDensities[i]

   #Load phi from GENE here for the specific run. Then split by mode.
   ampSpectra = readAmpSpectra(ampSpectraFiles[i])
   phi = np.array(ampSpectra[1][1])
   phi = phi[1:] #Drop krho = 0 term.
   rhostar_i = (vTi/cyc_i)/Lref
   phi = phi*rhostar_i*Ti/np.abs(z_e) #Normalize back to SI.

   #y = vPerp, x = vPar - will evaluate y integral first.
   for j, kthRho in enumerate(kthetaRhoi):
      if (j==1):
         break

      omega        = (omegaReal[j] + 1j*gamma[j])*(cyc_i/rhostar_i) #Convert from GENE to SI.
      oStarT_I     = lambda y, x: kthRho * (vT_I/c) * (rhoI/rho_i)  * ((-1/LnI)  - (1/2)*((y/vT_I)**2 + (x/vT_I)**2 - 3)*(1/L_TI)) #omega^T_star,I
      oStarT_i     = lambda y, x: kthRho * (vTi/c)  * (rho_i/rho_i) * ((-1/Ln_i) - (1/2)*((y/vTi)**2  + (x/vTi)**2  - 3)*(1/L_T))
      obarD_I      = lambda y, x: kthRho/R * (vT_I/c) * (rhoI/rho_i)  * ((1/2)*((y/vT_I)**2) + ((x/vT_I)**2)) * thetaIntegral #\bar{omega}_d,I
      obarDi       = lambda y, x: kthRho/R * (vTi/c)  * (rho_i/rho_i) * ((1/2)*((y/vTi)**2)  + ((x/vTi)**2))  * thetaIntegral
      omegaTerm    = lambda y, x: (np.conj(omega) - oStarT_I(y,x))/(np.conj(omega) - obarD_I(y,x)) + (omega - oStarT_I(y,x))/(omega - obarD_I(y,x))
      omegaTerm_i  = lambda y, x: (np.conj(omega) - oStarT_i(y,x))/(np.conj(omega) - obarDi(y,x))  + (omega - oStarT_i(y,x))/(omega - obarDi(y,x)) 
      xTerm_i      = lambda x:   scipy.exp(-1*(x/(2*vTi))**2)
      yTerm_i      = lambda y: y*scipy.exp(-1*(y/(2*vTi))**2)
      xTermI       = lambda x:   scipy.exp(-1*(x/(2*vT_I))**2)
      yTermI       = lambda y: y*scipy.exp(-1*(y/(2*vT_I))**2)
      besselTerm   = lambda y: scipy.special.jv(0, kthRho * (y/vTi)**2 * (cyc_i/cycI))
      besselTerm_i = lambda y: scipy.special.jv(0, kthRho * (y/vTi)**2 * (cyc_i/cyc_i))
      constantTermI  = kthRho * (z_i*zI/(c*vTi*m_i*T_I)) * phi[j]**2 * (mI/(2*np.pi*T_I))**(3/2) * nI
      constantTerm_i = kthRho * (z_i*z_i/(c*vTi*m_i*Ti)) * phi[j]**2 * (m_i/(2*np.pi*Ti))**(3/2) * n_i
   
      func     = lambda y, x: constantTermI*xTermI(x)*yTermI(y)*besselTerm(y)*omegaTerm(y,x)
      result   = complex_quadrature(func, -scipy.inf, scipy.inf, 0, scipy.inf)
      func_i   = lambda y, x: constantTerm_i*xTerm_i(x)*yTerm_i(y)*besselTerm_i(y)*omegaTerm_i(y,x)
      result_i = complex_quadrature(func_i, -scipy.inf, scipy.inf, 0, scipy.inf)

      print(np.real(result[0]))
      print(np.real(result_i[0]))

dirName = 'GoerlerImpurities'
if not os.path.exists(dirName):
    os.makedirs(dirName)
else:
   shutil.rmtree(dirName)
   os.makedirs(dirName)

#np.savez_compressed(dirName + '/flux.npz',       flux)
#np.savez_compressed(dirName + '/species.npz',    speciesNames[:-2])
#np.savez_compressed(dirName + '/growthData.npz', linearGrowthData)