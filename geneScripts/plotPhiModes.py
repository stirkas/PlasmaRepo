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

def plotFlux(dataOne, dataTwo, caseNumber):
   fig = plt.figure(num=None, figsize=(6,6), dpi=100)
   padding   = 10
   bigText   = 20
   smallText = 14
   plt.rcParams.update({'font.size': bigText})

   ax1 = fig.add_subplot(211)
   ax2 = fig.add_subplot(212, sharex=ax1)

   fig.text(0.5, 0.01, "t / ($L_{ref}$/$c_s$)", ha='center', fontsize=bigText)
   fig.text(0.01, 0.5, "|$\\phi$|", va='center', rotation='vertical', fontsize=bigText)

   t = np.zeros(np.shape(dataOne)[0])
   phiRealOne = np.zeros(np.shape(dataOne)[0])
   phiImagOne = np.zeros(np.shape(dataOne)[0])
   phiRealTwo = np.zeros(np.shape(dataTwo)[0])
   phiImagTwo = np.zeros(np.shape(dataTwo)[0])
   i = 0

   for datum in dataOne:
      t[i] = (np.abs(float(datum[0])))
      phiRealOne[i] = (np.abs(float(datum[1])))
      phiImagOne[i] = (np.abs(float(datum[2])))
      i = i + 1
   i = 0
   for datum in dataTwo:
      phiRealTwo[i] = (np.abs(float(datum[1])))
      phiImagTwo[i] = (np.abs(float(datum[2])))
      i = i + 1

   phiOne = np.sqrt(phiRealOne*phiRealOne + phiImagOne*np.conj(phiImagOne))
   phiTwo = np.sqrt(phiRealTwo*phiRealTwo + phiImagTwo*np.conj(phiImagTwo))
   
   p1 = ax1.plot(t, phiOne)
   p2 = ax1.plot(t, phiTwo)
   ax1.xaxis.set_tick_params(labelsize=smallText)
   ax1.yaxis.set_tick_params(labelsize=smallText)
   ax1.set_xlim(0,np.max(t))
 
   ax2.plot(t, phiOne)
   ax2.plot(t, phiTwo)
   ax2.xaxis.set_tick_params(labelsize=smallText)
   ax2.yaxis.set_tick_params(labelsize=smallText)
   ax2.set_yscale('log')
   ax2.set_ylim(.001, 10)
   ax2.set_xlim(0,np.max(t))

   ax1.grid()
   ax2.grid()
   #Add shared legend and shift things for it.
   #line_labels = ['$k_x$ = 0,                       $k_y$ = 6.36 [3*$k_{y,min}$]', '$k_x$ = 4.955 [4*$k_{x,min}$], $k_y$ = 0']
   #line_labels = ['$k_x$ = 4.955 [4*$k_{x,min}$], $k_y$ = 0']
   #fig.legend([p1,p2],labels=line_labels,bbox_to_anchor=(1,0), bbox_transform=ax2.transAxes, fontsize=smallText, framealpha=1)
   #fig.legend([p1,p2],labels=line_labels,loc='upper center', fontsize=smallText, framealpha=1)

   plt.tight_layout()
   #plt.show()
   plt.savefig('./geneScripts/GoerlerZonal/timeTracePlot_' + str(caseNumber) + '.pdf')

data1 = []
data2 = []
caseNumber = '11'
readPhiModes('./geneScripts/GoerlerZonal/timetraceelectrons_00' + str(caseNumber) + '_0.dat', data1)
readPhiModes('./geneScripts/GoerlerZonal/timetraceelectrons_00' + str(caseNumber) + '_1.dat', data2)
plotFlux(data1, data2, caseNumber)
