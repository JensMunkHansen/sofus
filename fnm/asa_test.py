import numpy as np
import matplotlib.pyplot as plt

import addpaths
import swig_fnm as fnm

from asa import asa

def create_time_string(seconds):
  m, s = divmod(seconds, 60)
  h, m = divmod(m, 60)
  ms = np.round(s * 1000) % 1000
  timeString = "%d:%02d:%02d,%03d" % (h, m, s, (1000*ms))
  return timeString;

plt.ion()

# Set up the array
nelex = 32
neley = 1
width = 0.245e-3
height = 7e-3
kerf_x = 0.03e-3
kerf_y = 0

f0 = 1e6
soundspeed = 1500
lambd = soundspeed / f0

ndiv = 10
a = fnm.ApertureFloat(nelex,width,kerf_x,height)
a.f0 = f0
a.c  = soundspeed
a.nDivW = ndiv
a.nDivH = ndiv
a.nthreads = 4

# Set up the coordinate grid
focus_x = 0
focus_y = 0
focus_z = 20 * lambd

dx = lambd/8
dy = lambd/8
dz = lambd/8
#dz = 2*dz

nx = 71
ny = 57
nz = 321

x = (np.r_[0:nx] - (nx-1)/2.0)*dx
y = (np.r_[0:ny] - (ny-1)/2.0)*dy
z = np.r_[0:nz]*dz

# Determine where the source pressure will be calculated
z0 = lambd/4
iy = ny/2

a.focus = [focus_x, focus_y, focus_z]

# Calculate the pressure

xs,ys = np.meshgrid(x,y,indexing='ij')

pos = np.c_[xs.flatten(), ys.flatten(), z0*np.ones(nx*ny)].astype(np.float32)

p0 = a.CalcCwFast(pos)
p0 = p0.reshape((nx,ny))

plt.figure()
plt.imshow(np.abs(p0),aspect='auto')

N=256

from timeit import default_timer as timer

start = timer()
p_asa = asa(p0,dx,dy,z,f0=f0,N=N)
end = timer()
print(create_time_string(end-start))

plt.figure()
plt.imshow(np.abs(p_asa[:,:,iy]),aspect='auto')

xs,zs = np.meshgrid(x,z)
pos = np.c_[xs.flatten(), y[iy]*np.ones(nx*nz), zs.flatten()].astype(np.float32)

p1 = a.CalcCwFast(pos)
p1 = p1.reshape((nz,nx))
plt.figure()
plt.imshow(np.abs(p1),aspect='auto')


