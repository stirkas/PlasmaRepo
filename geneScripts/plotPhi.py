#!/usr/bin/env python

import argparse
import os
import sys
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

from plotFlux import readFlux
from plotFlux import getFluxData

currentlySaving = False #Define globally so plotting routine sees if it changes.

def readPhi(fileName, data, version=2):
   readingData    = False
   reachedPhiData = False
   t = x = y = 0
   f = open(fileName, 'r')

   if (version >= 2): #My old GENE files looked different.
      for line in f.readlines():
         if (line.strip() == '#  phi'):
            reachedPhiData = True
            if (readingData and line.strip() == "#  phi"): #Ignore comments at top of file. Afterwards, phi token delineates start of new data set.
               t = t + 1
               x = 0
         else:
            if (reachedPhiData):
               if (t == 0): readingData = True #Reached first line of data.
               if (len(line.split()) > 0 and line[0] != '#'): #Ignore empty and comment lines. Not sure why there are x-values with no corresponding y sometimes.
                  data[t,x,:] = line.split()
               x = x + 1
   else:
      for line in f.readlines():
         if (line[0] == '#'):
            if (readingData): #Ignore comments at top of file. After, comments delineate timesteps.
               t = t + 1
               x = 0
         else:
            if (t == 0): readingData = True
            if (len(line.split()) > 0): #Ignore empty lines. Not sure why there are x-values with no corresponding y sometimes.
               data[t,x,:] = line.split()
            x = x + 1
   f.close()
   readingData    = False
   reachedPhiData = False

def update_anim(it,fig,phit,phikt,nt,xp,yp,kxp,kyp,plot,fluxData,tRatio,plotFlux=False):
   fig.clf()   
   ax1 = fig.add_subplot(211)
   ax2 = fig.add_subplot(212)
   ax1.clear()
   ax2.clear()
   phi = np.transpose(phit[it,:,:])
   phik = np.transpose(phikt[it,:,:])

   if (plotFlux and fluxData):
      fluxIt = it*tRatio if (it*tRatio < len(fluxData[0]) - 1) else len(fluxData[0]) - 1 #Flux timesteps usually different. Stop at max index.
      fluxItLabel = '$t_i$: ' + str(fluxIt)
      ax1.plot(fluxData[0], fluxData[1], label='$e^-$')
      ax1.plot(fluxData[0][fluxIt], fluxData[1][fluxIt], 'bo', label=fluxItLabel)
      ax1.title.set_text('Particle Flux')
      ax1.set_xlabel('$t$')
      ax1.set_ylabel('$\\langle Q_{ES}\\rangle$')
      ax1.grid()
      ax1.legend()
   else:
      im1 = ax1.contourf(xp,yp,phi,20,cmap='jet') #20 for colormap resolution.
      ax1.grid()
      ax1.title.set_text("$\\phi$")
      ax1.set_xlabel("$x/\\rho_i$")
      ax1.set_ylabel("$y/\\rho_i$")
      fig.colorbar(im1, ax=ax1)
   im2 = ax2.contourf(kxp,kyp,phik,20,cmap='jet')
   ax2.grid()
   ax2.title.set_text("$\\phi_k$")
   ax2.set_xlabel("$k_x\\rho_i$")
   ax2.set_ylabel("$k_y\\rho_i$")
   fig.colorbar(im2, ax=ax2)
   if (True): #Add axis limits here manually.
      ax2.set_xlim(-20,20)
      ax2.set_ylim(-20,20)
   plt.tight_layout()

   if (currentlySaving):
      print("Saving frame: " + str(it+1) + "/" + str(nt) + ".")
   else:
      print("Displaying frame: " + str(it+1) + "/" + str(nt) + ".")

      if ((it+1)==nt): #Since it index starts at 0, but nt is a count so doesn't include 0.
         if (plot):
            plt.close(fig)

def plotPhi(nt,nx,nky,lx,ly,rDataPath,kDataPath,fluxfile=None,tRatio=10,version=2,plot=False,save=False,savefile=None):
   nx   = nx    + 1 #I think GENE uses 2 gridpoints per mode + origin(0).
   nkyp = 2*nky - 1 #Num grid points in ky space. Not sure why - 1.
   nkx  = nx    - 1 #Not sure why GENE data has 1 less point for nkx.
   ny   = nky*2 + 1 #I think GENE uses 2 gridpoints per mode + origin(0).

   phit  = np.zeros((nt,nx,ny))
   phikt = np.zeros((nt,nkx,nkyp))
   xp = np.linspace(-lx/2, lx/2, nx)
   yp = np.linspace(-ly/2, ly/2, ny)
   kxp = (2*np.pi*nkx/lx)*np.linspace(-1/2 + 1/nkx,1/2,nkx) #For some reason x needs to be offset by 1.
   #Need to double the range for y. Apparently GENE does the full range of modes positive, then reflects.
   kyp = (2*np.pi*nky/ly)*np.linspace(-1,1,nkyp)
   
   readPhi(rDataPath, phit,  version)
   readPhi(kDataPath, phikt, version)

   fluxData = []
   plotFlux = False
   if (fluxfile != None and os.path.exists(fluxfile) and os.path.getsize(fluxfile) > 0):
      fluxType = 1 #0=g, 1=q, 2=p.
      readFlux(fluxfile, fluxData)
      fluxData = getFluxData(fluxType, fluxData)
      plotFlux = True

   # Set up formatting for the movie files
   fig = plt.figure()
   anim=animation.FuncAnimation(fig,update_anim,frames=nt,fargs=(fig,phit,phikt,nt,xp,yp,kxp,kyp,plot,fluxData,tRatio,plotFlux),repeat=False)
   if (plot):
      plt.show()
   if (save):
      global currentlySaving
      currentlySaving = True
      Writer = animation.writers['ffmpeg'] #Requires ffmpeg package on linux.
      writer = Writer(fps=15, bitrate=-1, codec='h264')
      if (savefile):
         anim.save(savefile, writer=writer)

plotPhi(1879, 192, 16, 5.63558, 2.96377,
        './GoerlerZonalETG/adIonPhi.dat', './GoerlerZonalETG/adIonPhik.dat', plot=True)
        #, fluxfile='./GoerlerZonalETG/kinIonFlux.dat',
        #save=True, savefile='./GoerlerZonalETG/phiZonalETGwFluxKinIon.mp4')

if __name__ == "__main__":
   parser = argparse.ArgumentParser("Script for plotting/saving GENE ASCII phi/phik output.")
   parser.add_argument("nt",  help="Number of timesteps.",     type=int)
   parser.add_argument("nx",  help="Number of x grid points.", type=int)
   parser.add_argument("nky", help="Number of ky modes.",      type=int)
   parser.add_argument("lx",  help="Box size in x.",           type=float)
   parser.add_argument("ly",  help="Box size in y.",           type=float)
   parser.add_argument("realDataPath",      help="Path to x-y data file.")
   parser.add_argument("kSpaceDataPath",    help="Path to kspace data file.")
   parser.add_argument("-p", "--plot",      help="Show plot over time.",               action="store_true")
   parser.add_argument("-s", "--save",      help="Save output plot. (15fps,h264,mp4)", action="store_true")
   parser.add_argument("-f", "--savefile",  help="Save filename.")
   parser.add_argument("-q", "--fluxfile",  help="Optional path to flux data file.")
   parser.add_argument("-t", "--timeRatio", help="Timestep ratio b/w flux data and phi data.")
   parser.add_argument("-v", "--version",   help="GENE version. Supports data phi files from 1 and 2.")
   args = parser.parse_args()

   #Require save and savefile together.
   if (len([x for x in (args.save,args.savefile) if (x is not None) and (x is not False)]) == 1):
      parser.error('--save and --savefile must be given together')

   plotPhi(args.nt, args.nx, args.nky, args.lx, args.ly,
           args.realDataPath, args.kSpaceDataPath, fluxfile=args.fluxfile, tRatio=args.timeRatio,
           version=args.version, plot=args.plot, save=args.save, savefile=args.savefile)