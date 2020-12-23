#!/usr/bin/env python

import argparse
import sys
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

def readPhiModes(fileName, data):
   f = open(fileName, 'r')

   for line in f.readlines():
      if (line[0] != '#'): #Ignore comment lines.
         if (len(line.split()) > 0): #Ignore empty lines. Not sure why there are x-values with no corresponding y sometimes.
            print(line.split())
            data.append(line.split())
   f.close()

def plotFlux(caseNumber, dataArray, labels, kmins):
   fig = plt.figure(num=None, figsize=(6,6), dpi=100)
   bigText   = 20
   smallText = 14
   smallerText = 10
   plt.rcParams.update({'font.size': bigText})
   
   ax1 = fig.add_subplot(211)
   ax2 = fig.add_subplot(212, sharex=ax1)

   #title = 'k$_{x,min}\\rho_i$=' + kmins[0] + ', k$_{y,min}\\rho_i$=' + kmins[1]
   #plt.suptitle('k$_{x,min}\\rho_i$='+kmins[0] + ', k$_{y,min}\\rho_i$='+kmins[1], fontsize=smallText, y=.995)

   fig.text(0.5, 0.01, "t / ($L_{ref}$/$c_s$)", ha='center', fontsize=smallText)
   fig.text(0.01, 0.5, "|$\\phi$|", va='center', rotation='vertical', fontsize=smallText)

   #Plot each mode separately. (First index is number of data sets concat'd)
   axisList = []
   for i in range(np.shape(dataArray)[0]):
      t       = np.zeros(np.shape(dataArray[i,:,:])[0]) #Get the data array sizes.
      phiReal = np.zeros(np.shape(dataArray[i,:,:])[0]) #Plot should throw an error if they dont all match through all loops.
      phiImag = np.zeros(np.shape(dataArray[i,:,:])[0])
      
      #Grab data for each mode plot j.
      j = 0
      for datum in dataArray[i,:,:]:
         t[j] = (np.abs(float(datum[0])))
         phiReal[j] = (np.abs(float(datum[1])))
         phiImag[j] = (np.abs(float(datum[2])))
         j = j + 1
   
      phi = np.sqrt(phiReal*phiReal + phiImag*np.conj(phiImag))
   
      #Plot data and save axis info for legend.
      axisList.append(ax1.plot(t, phi))
      ax2.plot(t, phi)
      
   ax1.xaxis.set_tick_params(labelsize=smallText)
   ax1.yaxis.set_tick_params(labelsize=smallText)
   ax1.set_xlim(0,np.max(t))
   ax2.xaxis.set_tick_params(labelsize=smallText)
   ax2.yaxis.set_tick_params(labelsize=smallText)
   ax2.set_yscale('log')
   ax2.set_ylim(.001, 10)
   ax2.set_xlim(0,np.max(t))

   ax1.grid()
   ax2.grid()
   #Add shared legend and shift things for it if necessary.
   #fig.legend([p1,p2],labels=labels,bbox_to_anchor=(1,0), bbox_transform=ax2.transAxes, fontsize=smallerText, framealpha=1)
   fig.legend(axisList,labels=labels,loc='upper right', fontsize=smallerText, framealpha=1)

   plt.tight_layout()
   plt.savefig('./modes' + str(caseNumber) + '.pdf')
   plt.show()

data1 = []
data2 = []
data3 = []
data4 = []

#Zonal flow folder scans.
#caseNumber = '15'
#readPhiModes('.zonalScan/kx0ky1_00' + str(caseNumber) + '.dat', data1)
#readPhiModes('.zonalScan/kx2ky0_00' + str(caseNumber) + '.dat', data2)
#readPhiModes('.zonalScan/kx4ky0_00' + str(caseNumber) + '.dat', data3)
#readPhiModes('.zonalScan/kx5ky0_00' + str(caseNumber) + '.dat', data4)
#labelOne    = 'k$_{y,i}$ = 1'
#labelTwo    = 'k$_{x,i}$ = 2'
#labelThree  = 'k$_{x,i}$ = 4'
#labelFour   = 'k$_{x,i}$ = 5'
#kmins = ['1.24','21.2'] #[kxmin, kymin]

#Original zonal scans.
caseNumber = 11
kmins = ['1.24', '2.12']
labelArr = np.array(['k$_x\\rho_i$=4.955', 'k$_y\\rho_i$=6.36'])
readPhiModes('./timetraceelectrons_0011_0.dat', data1)
readPhiModes('./timetraceelectrons_0011_1.dat', data2)
dataArr = np.array([data1, data2])

plotFlux(caseNumber, dataArr, labelArr, kmins)
