import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#Init spatial grid.
lx = 2*np.pi
nx = 256
dx = lx/nx
ly = 2*np.pi
ny = 256
dy = ly/ny

#Init temporal grid.
tFinal = 10
nt = 250
dt = tFinal/nt

#Create initial potential.
x = np.arange(nx)*dx
y = np.arange(ny)*dy
X,Y = np.meshgrid(x,y)

s=2
s2=s**2
r=(X-lx/2)**2+(Y-ly/2)**2
phi = np.exp(-r/s2)

plt.imshow(phi, origin='lower', interpolation='none')
plt.colorbar()
plt.show()