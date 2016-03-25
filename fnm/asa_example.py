import sys
import numpy as np
from timeit import default_timer as timer
import matplotlib.pyplot as plt

import addpaths
import swig_fnm as fnm
from asa import asa

sys.stdout.write('=================================[ ASA_Example.py ]================================\n\n');
sys.stdout.write('This script calculates the 3D pressure of a 128-element array of rectangular\n');
sys.stdout.write('transducers by using the Angular Spectrum Approach to propagate an initial field\n');
sys.stdout.write('calculated with the Fast Nearfield Method. The script produces two plots:\n');
sys.stdout.write('the initial FNM pressure and the pressure at y = 0 as calculated with the ASA.\n\n');

plt.ion()

ele_x = 128
ele_y = 1
width = 0.245e-3
height = 7e-3
kerf_x = 0.03e-3
kerf_y = 0
d = ele_x * (width+kerf_x) # Array aperture
N = 2 # f-number
# Create planar array of rectangular transducers

xdcr_array = fnm.ApertureFloat(ele_x, width, kerf_x, height)
# Use lossless medium
#medium = set_medium('lossless')
c = 1500
# Center frequency and wavelength
f0 = 3e6
xdcr_array.f0 = f0

lambd = c/f0

# Set the focus to be at the desired f-number
focus_x = 0
focus_y = 0
focus_z = d * N

# Set up the coordinate grid
xmin = -((ele_x/2.0) * (width+kerf_x))*1.5
xmax = -xmin
ymin = -((ele_y/2.0) * (height+kerf_y))*1.5
ymax = -ymin
zmin = lambd/4
zmax = focus_z*2

dx = lambd/2
dy = lambd/2
dz = lambd/2

x = np.r_[xmin:xmax+0.5*dx:dx]
y = np.r_[ymin:ymax+0.5*dy:dy]
z = np.r_[zmin:zmax+0.5*dz:dz]

# Determine where the source pressure will be calculated
z0 = lambd/4
y_index = np.floor((ymax-ymin)/2/dy)
# Coordinate grids to calclate the initial pressure (x-y plane) and final
# pressure (x-z plane)
#cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0)
#cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax)

# Focus the array
ndiv = 200
xdcr_array.nDivH = ndiv
xdcr_array.nDivW = ndiv
xdcr_array.focus = [focus_x,focus_y,focus_z]

# Calculate the pressure
ndiv = 10
xdcr_array.nDivH = ndiv
xdcr_array.nDivW = ndiv
sys.stdout.write('Calculating initial pressure plane with FNM... ')
start = timer()

[xs,ys] = np.meshgrid(x,y,indexing='ij')
pos = np.c_[xs.flatten(), ys.flatten(), z0*np.ones(xs.shape).flatten()].astype(np.float32)

#p0 = cw_pressure(xdcr_array,cg_p0,medium,ndiv,f0,'fnm sse')
p0 = xdcr_array.CalcCwFast(pos)
p0 = p0.reshape((len(x),len(y)))
p0 = p0.astype(np.complex64)
sys.stdout.write('done in %f s.\n' % (timer() - start))

sys.stdout.write('Calculating 3D pressure (%d points) with ASA... ' % (len(x) * len(y) * len(z)))
start = timer()
Nx=len(x)
Ny=len(y)
Nx=1024
Ny=1024
#p_asa = cw_angular_spectrum(p0,cg_3d,medium,f0,1024,'Pa')
p_asa = asa(p0,dx,dy,z,Nx=Nx,Ny=Ny,f0=f0,c=c,beta=0.0)
sys.stdout.write('done in %f s.\n' % (timer() - start))

# Show the initial pressure
plt.figure(1)
plt.imshow(np.abs(p0).T,extent=1000.0*np.r_[x[0],x[-1],y[0],y[-1]],aspect='auto',interpolation='nearest')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.title("p0 (Calculated with FNM at z = %3.2f mm)" % (z0*1000))
plt.show()

# Show the 3D field calculated with ASA
plt.figure(2)
plt.imshow(np.abs(np.squeeze(p_asa[:,:,y_index])),extent=1000.0*np.r_[x[0],x[-1],z[-1],z[0]],aspect='auto',interpolation='nearest')
plt.xlabel('x (mm)')
plt.ylabel('z (mm)')
plt.title('ASA Pressure (y=0)')
plt.show()
