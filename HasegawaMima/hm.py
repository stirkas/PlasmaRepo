import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

#Init spatial vars.
lx = 10
nx = 64
dx = lx/nx
ly = 10
ny = 64
mx = my = 8
dy = ly/ny
kappa = .1

#Init temporal grid.
nt = 400
dt = 2e-2 #Timestep taken from Numata code.

#Create initial grids.
x = np.arange(nx)*dx
y = np.arange(ny)*dy
X,Y = np.meshgrid(x,y)

#Create grid vals for fourier space. Note, drop end of linspace to match the fft of the arange in real space.
kx = nx*np.linspace(-1/2,1/2,nx,endpoint=False)
ky = ny*np.linspace(-1/2,1/2,ny,endpoint=False)
#Shift to match fft output for calculations.
kx = np.fft.ifftshift(kx)
ky = np.fft.ifftshift(ky)
KX,KY = np.meshgrid(kx,ky)
#Gather up alias vectors. They remove outer 1/3 of k-modes (on each side) to keep nonlinear vals in our k-range (2/3 rule)
#Store in nonshifted mode to use for calculations.
kxd = np.r_[np.ones(nx//3),np.zeros(nx//3+nx%3),np.ones(nx//3)]
kyd = np.r_[np.ones(ny//3),np.zeros(ny//3+ny%3),np.ones(ny//3)]
KXD,KYD = np.meshgrid(kxd,kyd)

#Create a phi and take its fft. Shifting values in order of negative to positive.
#waveFreq = 16
#phi = np.cos(waveFreq*2*np.pi*Y/ly)*np.cos(waveFreq*2*np.pi*X/lx)
#phi = np.exp(-((X-lx/2)**2 + (Y-ly/2)**2)/4)
waveFreq = 5
A=np.zeros((mx,my))+0.*1j
phi=np.zeros((nx,ny))+0.*1j

A[0,2]=1.*np.exp(2.1J)
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

phik = np.fft.fft2(phi)

#Allocate space for time routine.
saveRate = 4
numFrames = nt//saveRate #Save images once every twenty five dt's
phit  = np.zeros((numFrames,nx,ny))
phikt = np.zeros((numFrames,nx,ny))
phit[0,:,:] = np.real(phi)
phikt[0,:,:] = np.abs(np.fft.fftshift(phik))

#Set up kx^2
kx2 = kx**2
ky2 = ky**2
KX2, KY2 = np.meshgrid(kx2,ky2)
kconst = 1/(1+KX2+KY2) #Save off const for later calculations.

def adv(phik):
   #zeta = del2(phi)
   zetak = -(KX2 + KY2)*phik

   #Define derivs for main equation.
   phikx  = 1j*KX*phik;   phix = np.real(np.fft.ifft2(phikx*KXD*KYD))
   phiky  = 1j*KY*phik;   phiy = np.real(np.fft.ifft2(phiky*KXD*KYD))
   zetakx = 1j*KX*zetak; zetax = np.real(np.fft.ifft2(zetakx*KXD*KYD))
   zetaky = 1j*KY*zetak; zetay = np.real(np.fft.ifft2(zetaky*KXD*KYD))

   derivative = kconst*(np.fft.fft2(zetax*phiy - zetay*phix) + kappa*phiky)
   
   return derivative

def update_anim(it):
   fig.clf()
   ax1 = fig.add_subplot(121)
   ax2 = fig.add_subplot(122)
   ax1.clear()
   ax2.clear()
   im1 = ax1.contourf(X,Y,phit[it,:,:])
   im2 = ax2.contourf(np.fft.fftshift(KX), np.fft.fftshift(KY), phikt[it,:,:])
   ax1.grid()
   ax2.grid()
   plotRatio = 3/2
   ax2.set_xlim(-plotRatio*waveFreq, plotRatio*waveFreq)
   ax2.set_ylim(-plotRatio*waveFreq, plotRatio*waveFreq)
   fig.colorbar(im1, ax=ax1)
   fig.colorbar(im2, ax=ax2)
   plt.tight_layout()

#Main loop
for it in range(1,nt):
   rk1 = adv(phik)
   rk2 = adv(phik + 0.5*dt*rk1)
   rk3 = adv(phik + 0.5*dt*rk2)
   rk4 = adv(phik + dt*rk3)

   rk = (rk1 + 2*rk2 + 2*rk3 + rk4)/6

   phik = phik + dt*rk #RK4 advance.
   
   #Store plots as specified.
   if ((it%(saveRate))==0):
      print("Data slot: " + str(it))
      phit[it//saveRate,:,:]  = np.real(np.fft.ifft2(phik))
      phikt[it//saveRate,:,:] = np.abs(np.fft.fftshift(phik))

fig = plt.figure()
anim=animation.FuncAnimation(fig,update_anim,frames=numFrames,repeat=False)
plt.show()