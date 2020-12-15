#!/usr/bin/env python

import argparse
import sys
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

def readFlux(fileName, data):
   f = open(fileName, 'r')

   for line in f.readlines():
      if (line[0] != '#'): #Ignore comment lines.
         if (len(line.split()) > 0): #Ignore empty lines. Not sure why there are x-values with no corresponding y sometimes.
            data.append(line.split())
   f.close()

def getFluxData(fluxVal, data):
   t = [] #time
   g = [] #particle flux
   q = [] #heat flux
   p = [] #momentum flux
   outputData = []

   #Set up data arrays.
   for dataArray in data:
      t.append(float(dataArray[0]))
      g.append(float(dataArray[1]))
      q.append(float(dataArray[2]))
      p.append(float(dataArray[3]))

   if (fluxVal == 0):
      outputData = g
   elif (fluxVal == 1):
      outputData = q
   elif (fluxVal == 2):
      outputData = p

   return [t, outputData]

def plotFlux(dataPath,fluxVal,title='FluxData',plot=False,save=False,savefile=None):
   data = []
   readFlux(dataPath, data)
   [t, outputData] = getFluxData(fluxVal, data)
   
   plt.plot(t, outputData)
   
   plt.title(title, fontsize=16)
   plt.xlabel('$t\\;/\\;(L_{ref}/c_s)$', fontsize=16)
   plt.ylabel('$\\langle Q_{ES}\\rangle$', fontsize=16)
   plt.grid()
   plt.tight_layout()
   
   if (plot):
      plt.show()
   if (save):
      plt.savefig(savefile + '.pdf')

def getFluxAverage(data):
   average = 0
   for datum in data:
      average = average + datum
   average = average/len(data)

   return average

def getFluxStdDev(data):
   stdDev = 0
   avg = getFluxAverage(data)

   for datum in data:
      stdDev = stdDev + (datum - avg)**2
   stdDev = stdDev/len(data)
   stdDev = np.sqrt(stdDev)

   return stdDev

if __name__ == "__main__":
   parser = argparse.ArgumentParser("Script for plotting/saving GENE ASCII flux output. Add x,y labels manually.")
   parser.add_argument("dataPath",         help="Path to data file.")
   fv = parser.add_argument("fluxVal",     help="0=particle, 1=heat, 2=momentum", type=int)
   parser.add_argument("-t", "--title",    help="Plot title. Involving particle species for instance.")
   parser.add_argument("-p", "--plot",     help="Show plot over time.",                      action="store_true")
   parser.add_argument("-s", "--save",     help="Save output file (pdf for lossless image)", action="store_true")
   parser.add_argument("-f", "--savefile", help="Save filename.", type=str)
   args = parser.parse_args()

   #Require save and savefile together.
   if (len([x for x in (args.save,args.savefile) if (x is not None) and (x is not False)]) == 1):
      parser.error('--save and --savefile must be given together')

   if (args.fluxVal < 0 or args.fluxVal > 2):
      parser.error('Unsupported flux type: ' + str(args.fluxVal) + '. Please try one of the following: ' +
       fv.help + '.')

   plotFlux(args.dataPath, args.fluxVal, title=args.title, plot=args.plot, save=args.save, savefile=args.savefile)