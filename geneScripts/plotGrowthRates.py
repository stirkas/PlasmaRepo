#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

def readGrowthRates(fileName):
   data = []
   f = open(fileName, 'r')

   for line in f.readlines():
      if (line[0] != '#'): #Ignore comment lines.
         if (len(line.split()) > 0): #Ignore empty lines. Not sure why there are x-values with no corresponding y sometimes.
            line = line.split() #Remove index from front as well.
            line[:] = [x for x in line if x != "|"] #Remove all dividers from file.
            data.append(line[1:])
   f.close()

   kyRhoi = []
   gamma  = []
   omega  = []

   #Set up data arrays.
   for dataArray in data:
      kyRhoi.append(float(dataArray[0]))
      gamma.append(float(dataArray[1]))
      omega.append(float(dataArray[2]))
      #TODO: Add more modes if necessary.

   return [kyRhoi, gamma, omega]

if __name__ == "__main__":
   runs = ['adElPure', 'adElProt', 'adElNeon', 'kinElPure', 'kinElProt', 'kinElTung']
   dataFiles = ['./GoerlerImpLin/scanfiles0000/scan.log', './GoerlerImpLin/scanfiles0001/scan.log',
                './GoerlerImpLin/scanfiles0002/scan.log', './GoerlerImpLin/scanfiles0004/scan.log',
                './GoerlerImpLin/scanfiles0005/scan.log', './GoerlerImpLin/scanfiles0007/scan.log']

   runs = runs[3:]
   dataFiles = dataFiles[3:]

   fig,axs = plt.subplots(2,1)

   for i, run in enumerate(runs):
      data = []
      [kyRhoi, gamma, omega] = readGrowthRates(dataFiles[i])
   
      scaleFactor = 1.05
      axs[0].plot(kyRhoi, gamma, marker='*', label=runs[i])
      axs[0].set_xlabel('k$_y$$\\rho_i$')
      axs[0].set_ylabel('$\\gamma$')
      axs[0].grid()
      axs[1].plot(kyRhoi, omega, marker='*', label=runs[i])
      axs[1].set_xlabel('k$_y$$\\rho_i$')
      axs[1].set_ylabel('$\\omega$')
      axs[1].grid()
   
   plt.legend(loc='upper right')
   plt.tight_layout()
   plt.show()