import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import random
import time

print("Initializing variables.")

#Init spatial vars.
mRat = (1/2000) #m_e/m_i
lx = 5.63558 #(2*np.pi/.15)*np.sqrt(mRat)
nx = 193
dx = lx/nx
ly = 2.96377 #(2*np.pi/.15)*np.sqrt(mRat)
ny = 33
mx = my = 8
dy = ly/ny

tau = 1 #T_e/T_i
eta = 3 #r_n/r_t #Usually bigger than 1.
rnByRhoE = 500*np.sqrt(1/mRat) #r_n/rho_e, 500 = Haotian's r_n/rho_i

#Init temporal grid.
nt = 100000
dt = 0.0046332372 #(1/rnByRhoE) #rho_e/r_n

#Create initial grids.
x = np.arange(nx)*dx
y = np.arange(ny)*dy
X,Y = np.meshgrid(x,y)

#Create grid vals for fourier space. Note, drop end of linspace to match the fft of the arange in real space.
kx = (2*np.pi*nx/lx)*np.linspace(-1/2 + (1/nx),1/2,nx,endpoint=False)
ky = (2*np.pi*ny/ly)*np.linspace(-1/2 + (1/ny),1/2,ny,endpoint=False)
#Shift to match fft output for calculations.
kx = np.fft.ifftshift(kx)
ky = np.fft.ifftshift(ky)

KX,KY = np.meshgrid(kx,ky)
#Gather up alias vectors. They remove outer 1/3 of k-modes (on each side) to keep nonlinear vals in our k-range (2/3 rule)
#Store in nonshifted mode to use for calculations.
kxd = np.r_[np.ones(nx//3),np.zeros(nx//3+nx%3),np.ones(nx//3)]
kyd = np.r_[np.ones(ny//3),np.zeros(ny//3+ny%3),np.ones(ny//3)]
KXD,KYD = np.meshgrid(kxd,kyd)

#Set up vars specific to routines.
phi = np.zeros((nx,ny))
plotRatio = 1
waveFreq = 1
initialCase = 5
showPlot = False
saveAnim = True
currentlySaving = False #Turn on once files are being saved.
caseString = "Unspecified"
#Allocate initial conditions.
if (initialCase == 1):
   #2 strong modes.
   caseString = "TwoStrongModes"
   waveFreq = 16
   plotRatio = 12
   phi = np.cos(waveFreq*2*np.pi*Y/ly)*np.cos(waveFreq*2*np.pi*X/lx)
elif (initialCase == 2):
   #Gaussian
   caseString = "Gaussian"
   waveFreq = 5
   plotRatio = 3/2
   phi = np.exp(-((X-lx/2)**2 + (Y-ly/2)**2)/4)
elif (initialCase == 3):
   #Random strong mode + random weaker modes + random phase shifts in each.
   caseString = "StrongModeWithWeakerPerturbations"
   waveFreq = .75
   plotRatio = 3/2 * np.sqrt(1/mRat)
   A=np.zeros((mx,my))+0.*1j
   phi=np.zeros((nx,ny))+0.*1j
   A[0,2]=1*np.exp(2.1J)
   A[0,3]=0.1*np.exp(1.7J)
   A[7,3]=0.05*np.exp(1.1J)
   A[5,2]=0.05*np.exp(0.1J)
   A[3,6]=0.05*np.exp(0.7J)
   A[4,7]=0.05*np.exp(2.7J)
   for i in range(nx):
       for j in range(ny):
           for m1 in range(mx):
               for m2 in range(my):
                   phi[i,j]= phi[i,j]+A[m1,m2]*np.exp(1j*kx[m1]*x[i]+1j*ky[m2]*y[j])
   phi=np.transpose(np.real(phi))
   scaleFactor = 1 #.6
   phi = scaleFactor*phi
elif (initialCase == 4): #Fredys ICs
   mx = 16
   waveFreq = .4
   plotRatio = 5
   px=np.linspace(-mx/2,mx/2,mx)*2.*np.pi/lx
   px[8]=0
   px[7]=0
   py=np.arange(my)*2.*np.pi/ly
   A=np.zeros((mx,my))+0.*1j
   phi=np.zeros((nx,ny))+0.*1j
   # had at 0.3-0.5 loop for 64 nx.
   # had 0-8 for both modes

   for k in range(mx):
       rand=random.uniform(0,2*np.pi)
       A[int(random.uniform(4,12)),int(random.uniform(0,8))]=random.uniform(0.08,0.15)*np.exp(1j*rand)

   A[8,int(random.uniform(3*my/8,4*my/8))] = 0.11*np.exp(1j*random.uniform(0,2*np.pi))

   rand=random.uniform(0,2*np.pi)
   A[8,int(my/4)]=1.*np.exp(1j*rand)  


   # Actual IC

   for i in range(nx):
       for j in range(ny):
           for m1 in range(mx):
               for m2 in range(my):
                   phi[i,j]=phi[i,j]+A[m1,m2]*np.exp(1j*px[m1]*x[i]+1j*py[m2]*y[j])
   phi=np.real(phi)
   phi=np.transpose(phi)
elif (initialCase == 5):
   #Load output from GENE.
   #Begin loading data.
   fileName = '/home/stirkas/Workspace/GENE/phi_0021_0025_r.dat'
   readingData = False
   f = open(fileName, 'r')
   t = x = y = 0
   phit = np.zeros((637, 33, 193)) #[t,x,y]

   for line in f.readlines():
      if (line[0] == '#'):
         if (readingData): #Ignore comments at top of file. After, comments delineate timesteps.
            t = t + 1
            x = 0
      else:
         if (t == 0): readingData = True
         if (len(line.split()) > 0): #Ignore empty lines. Not sure why there are x-values with no corresponding y sometimes.
            vals = line.split()
            phit[t,:,x] = vals
         x = x+1

   f.close()

   phi = phit[470,:,:]*(10**-2)

else:
   sys.exit("Invalid initial conditions. Current value: " + caseString + ".")

#Normalize some things nicely for ETG.
phi = phi*10**-4
phik = np.fft.fft2(phi)
print("Loaded initial conditions: " + caseString)

#Allocate space for time routine.
saveRate = 50
print("Running for " + str(nt) + " frames and saving data every " + str(saveRate) + " frames.")
numFrames = nt//saveRate #Save images once every twenty five dt's
phit  = np.zeros((numFrames,ny,nx))
phikt = np.zeros((numFrames,ny,nx))
phit[0,:,:]  = np.real(phi)
phikt[0,:,:] = np.abs(np.fft.fftshift(phik))

#Set up kx^2
kx2 = kx**2
ky2 = ky**2
KX2, KY2 = np.meshgrid(kx2,ky2)
kconst = 1/(1-(1+tau)*(KX2+KY2)/(2*tau)) #Save off const for later calculations.

def adv(phik):
   zetak = -(KX2+KY2)*phik

   #Define derivs for main equation.
   phikx  = 1j*KX*phik;   phix = np.real(np.fft.ifft2(phikx*KXD*KYD))
   phiky  = 1j*KY*phik;   phiy = np.real(np.fft.ifft2(phiky*KXD*KYD))
   zetakx = 1j*KX*zetak; zetax = np.real(np.fft.ifft2(zetakx*KXD*KYD))
   zetaky = 1j*KY*zetak; zetay = np.real(np.fft.ifft2(zetaky*KXD*KYD))

   term1 = ((1+tau)**2)*rnByRhoE/(4*(tau**2))
   term2 = (1+eta)/(2*tau)
   term3 = (1+tau)*(1+eta)/(4*tau)
   derivative = kconst*(term1*np.fft.fft2(phix*zetay-zetax*phiy) + term2*phiky + term3*zetaky)

   #print('---')
   #print(np.amax(phik))
   #print(np.amax(phikx))
   #print(np.amax(zetak))
   #print(np.amax(zetakx))
   #print(np.amax(derivative))
   #time.sleep(.1)

   return derivative

#Main loop
print("Starting main data loop.")
for it in range(1,nt+1):
   rk1 = adv(phik)
   rk2 = adv(phik + 0.5*dt*rk1)
   rk3 = adv(phik + 0.5*dt*rk2)
   rk4 = adv(phik + dt*rk3)

   rk = (rk1 + 2*rk2 + 2*rk3 + rk4)/6

   phik = phik + dt*rk #RK4 advance.
   
   #Store plots as specified.
   if ((it%(saveRate))==0):
      print("Storing frame data: " + str(it) + "/" + str(nt) + ".")
      phit[(it//saveRate)-1,:,:]  = np.real(np.fft.ifft2(phik))
      phikt[(it//saveRate)-1,:,:] = np.abs(np.fft.fftshift(phik))

print("Finished storing run data.")

def update_anim(it):
   fig.clf()   
   ax1 = fig.add_subplot(121)
   ax2 = fig.add_subplot(122)
   ax1.clear()
   ax2.clear()
   im1 = ax1.contourf(X,Y, phit[it,:,:], 20, cmap='jet')
   im2 = ax2.contourf(np.fft.fftshift(KX), np.fft.fftshift(KY), phikt[it,:,:], 20, cmap='jet')
   ax1.grid()
   ax2.grid()
   ax1.title.set_text("$\\phi$")
   ax2.title.set_text("$\\phi_k$")
   plotSize = 20
   ax2.set_xlim(-plotSize, plotSize)
   #ax2.set_xlim(-plotRatio*waveFreq, plotRatio*waveFreq)
   ax1.set_xlabel("x")
   ax2.set_xlabel("$k_x$")
   ax2.set_ylim(-plotSize, plotSize)
   #ax2.set_ylim(-plotRatio*waveFreq, plotRatio*waveFreq)
   ax1.set_ylabel("y")
   ax2.set_ylabel("$k_y$")
   fig.colorbar(im1, ax=ax1)
   fig.colorbar(im2, ax=ax2)
   plt.tight_layout()

   if (currentlySaving):
      print("Saving frame: " + str(it+1) + "/" + str(numFrames) + ".")
   else:
      print("Displaying frame: " + str(it+1) + "/" + str(numFrames) + ".")

   if ((it+1)==numFrames): #Since it index starts at 0, but numFrames is a count so doesn't include 0.
      if (showPlot):
         plt.close(fig)

print("Starting animation.")

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg'] #Requires ffmpeg package on linux.
writer = Writer(fps=30, bitrate=-1, codec='h264')

fig = plt.figure(num=None, figsize=(10,5))
anim=animation.FuncAnimation(fig,update_anim,frames=numFrames,repeat=False)

if (showPlot):
   plt.show()

if (saveAnim):
   print("Saving animation. This will probably take a few minutes...")
   currentlySaving = True
   anim.save('HasegawaMima/hmETG_' + str(initialCase) + '.mp4', writer=writer)
sys.exit("Animation complete.")