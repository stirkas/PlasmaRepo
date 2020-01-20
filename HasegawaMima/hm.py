import numpy as np
import matplotlib.pyplot as plt

#Init spatial grid.
lx = 10
nx = 256
dx = lx/nx
ly = 10
ny = 256
dy = ly/ny

#Init temporal grid.
tFinal = 10
nt = 250
dt = tFinal/nt

#Create initial grids.
x = np.arange(nx)*dx
y = np.arange(ny)*dy
X,Y = np.meshgrid(x,y)

#Create a phi and take its fft. Shifting values in order of negative to positive.
phi = np.cos(10*2*np.pi*Y/ly) + np.cos(4*2*np.pi*X/lx)
phik = np.fft.fft2(phi)/(nx*ny)
phik = np.fft.fftshift(phik)

#Create grid vals for fourier space. Note, drop end of linspace to match the fft of the arange in real space.
phifx = nx*np.linspace(-1/2,1/2,nx,endpoint=False)
phify = ny*np.linspace(-1/2,1/2,ny,endpoint=False)

fig, axs = plt.subplots(1,2)
im0 = axs[0].contourf(X, Y, phi)
im1 = axs[1].contourf(phifx, phify, np.abs(phik))
fig.colorbar(im0, ax=axs[0])
fig.colorbar(im1, ax=axs[1])
axs[0].grid()
axs[1].grid()
plt.tight_layout()
plt.show()