#!/usr/bin/env python

import argparse
import sys
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

currentlySaving = False #Define globally so plotting routine sees if it changes.

def readPhi(fileName, data):
   readingData = False
   t = x = y = 0
   f = open(fileName, 'r')

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
   readingData = False

def update_anim(it,fig,phit,phikt,nt,xp,yp,kxp,kyp,plot):
   fig.clf()   
   ax1 = fig.add_subplot(211)
   ax2 = fig.add_subplot(212)
   ax1.clear()
   ax2.clear()
   phi = np.transpose(phit[it,:,:])
   phik = np.transpose(phikt[it,:,:])
   im1 = ax1.contourf(xp,yp,phi,20,cmap='jet') #20 for colormap resolution.
   im2 = ax2.contourf(kxp,kyp,phik,20,cmap='jet')
   ax1.grid()
   ax2.grid()
   ax1.title.set_text("$\\phi$")
   ax2.title.set_text("$\\phi_k$")
   ax1.set_xlabel("$x/\\rho_i$")
   ax2.set_xlabel("$k_x\\rho_i$")
   ax1.set_ylabel("$y/\\rho_i$")
   ax2.set_ylabel("$k_y\\rho_i$")
   fig.colorbar(im1, ax=ax1)
   fig.colorbar(im2, ax=ax2)
   plt.tight_layout()

   if (currentlySaving):
      print("Saving frame: " + str(it+1) + "/" + str(nt) + ".")
   else:
      print("Displaying frame: " + str(it+1) + "/" + str(nt) + ".")

      if ((it+1)==nt): #Since it index starts at 0, but nt is a count so doesn't include 0.
         if (plot):
            plt.close(fig)

def plotPhi(nt,nx,nky,lx,ly,rDataPath,kDataPath,plot=False,save=False,savefile=None):
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
   
   readPhi(rDataPath, phit)
   readPhi(kDataPath, phikt)

   # Set up formatting for the movie files
   fig = plt.figure()
   anim=animation.FuncAnimation(fig,update_anim,frames=nt,fargs=(fig,phit,phikt,nt,xp,yp,kxp,kyp,plot),repeat=False)
   if (plot):
      plt.show()
   if (save):
      global currentlySaving
      currentlySaving = True
      Writer = animation.writers['ffmpeg'] #Requires ffmpeg package on linux.
      writer = Writer(fps=15, bitrate=-1, codec='h264')
      if (savefile):
         anim.save(savefile, writer=writer)

if __name__ == "__main__":
   parser = argparse.ArgumentParser("Script for plotting/saving GENE ASCII phi/phik output.")
   parser.add_argument("nt",  help="Number of timesteps.",     type=int)
   parser.add_argument("nx",  help="Number of x grid points.", type=int)
   parser.add_argument("nky", help="Number of ky modes.",      type=int)
   parser.add_argument("lx",  help="Box size in x.",           type=float)
   parser.add_argument("ly",  help="Box size in y.",           type=float)
   parser.add_argument("realDataPath",     help="Path to x-y data file.")
   parser.add_argument("kSpaceDataPath",   help="Path to kspace data file.")
   parser.add_argument("-p", "--plot",     help="Show plot over time.",               action="store_true")
   parser.add_argument("-s", "--save",     help="Save output plot. (15fps,h264,mp4)", action="store_true")
   parser.add_argument("-f", "--savefile", help="Save filename.")
   args = parser.parse_args()

   #Require save and savefile together.
   if (len([x for x in (args.save,args.savefile) if (x is not None) and (x is not False)]) == 1):
      parser.error('--save and --savefile must be given together')

   plotPhi(args.nt, args.nx, args.nky, args.lx, args.ly,
      args.realDataPath, args.kSpaceDataPath, plot=args.plot, save=args.save, savefile=args.savefile)