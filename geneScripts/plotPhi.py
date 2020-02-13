import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

nt = 509
nx = 194
ny = 33
nkx = nx-2 #GENE seems to drop size by 2 for k-space. 
nky = ny-2 #GENE seems to drop size by 2 for k-space.
lx = 5.63558
ly = 2.96377
t = x = y = 0
phit = np.zeros((nt,nx,ny)) #[t,x,y]
phikt = np.zeros((nt,nkx,nky))
readingData = False
currentlySaving = False
showPlot = False
xp = np.linspace(-lx/2, lx/2, nx)
yp = np.linspace(-ly/2, ly/2, ny)
xkp = (2*np.pi*nkx/lx)*np.linspace(-1/2+(1/nkx),1/2,nkx) #For some reason x needs to be offset by 1.
ykp = (2*np.pi*nky/ly)*np.linspace(-1/2,  1/2,nky)
#Begin loading data.
fileName = '/home/stirkas/Workspace/GENE/Output/phi_21-24_real.dat'
f = open(fileName, 'r')

for line in f.readlines():
   if (line[0] == '#'):
      if (readingData): #Ignore comments at top of file. After, comments delineate timesteps.
         t = t + 1
         x = 0
   else:
      if (t == 0): readingData = True
      if (len(line.split()) > 0): #Ignore empty lines. Not sure why there are x-values with no corresponding y sometimes.
         phit[t,x,:] = line.split()
      x = x + 1
f.close()

readingData = False
t = x = y = 0
fileName = '/home/stirkas/Workspace/GENE/Output/phi_21-24_k.dat'
f = open(fileName, 'r')

for line in f.readlines():
   if (line[0] == '#'):
      if (readingData): #Ignore comments at top of file. After, comments delineate timesteps.
         t = t + 1
         x = 0
   else:
      if (t == 0): readingData = True
      if (len(line.split()) > 0): #Ignore empty lines. Not sure why there are x-values with no corresponding y sometimes.
         phikt[t,x,:] = line.split()
      x = x + 1
f.close()

def update_anim(it):
   fig.clf()   
   ax1 = fig.add_subplot(211)
   ax2 = fig.add_subplot(212)
   ax1.clear()
   ax2.clear()
   phi = np.transpose(phit[it,:,:])
   phik = np.transpose(phikt[it,:,:])
   im1 = ax1.contourf(xp,yp,phi,20,cmap='jet') #20 for colormap resolution.
   im2 = ax2.contourf(xkp,ykp,phik)
   ax1.grid()
   ax2.grid()
   ax1.title.set_text("$\\phi$")
   ax2.title.set_text("$\\phi_k$")
   ax2.set_xlim(-40, 40)
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
      if (showPlot):
         plt.close(fig)

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg'] #Requires ffmpeg package on linux.
writer = Writer(fps=15, bitrate=-1, codec='h264')

fig = plt.figure()
anim=animation.FuncAnimation(fig,update_anim,frames=nt,repeat=False)

if (showPlot):
   plt.show()

currentlySaving = True
anim.save('geneScripts/phi.mp4', writer=writer)
