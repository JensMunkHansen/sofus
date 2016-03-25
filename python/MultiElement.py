import numpy as np

from fnm import (rect,linear_array)
from timeit import default_timer as timer

import matplotlib.pyplot as plt
plt.ion()

def create_time_string(seconds):
  m, s = divmod(seconds, 60)
  h, m = divmod(m, 60)
  ms = np.round(s * 1000) % 1000
  timeString = "%d:%02d:%02d,%03d" % (h, m, s, (1000*ms))
  return timeString;

print('This script uses the Fast Nearfield Method to calculate the CW pressure field of');
print('an array of 10 rectangular elements focused at a single point. The script');
print('outputs the pressure field.\n');

f0 = 1e6 # excitation frequency,Hz
soundspeed = 1500 # m/s
lamda = soundspeed / f0 # wavelength, m


#define a transducer structure/array
nelex = 10
neley = 1
kerf = 5.0e-4

width = 3e-3 # transducer width, m
height = 50e-3 # transducer height, m

d = nelex * (width+kerf)

xmin = -1.5 * d/2
xmax = 1.5 * d/2
ymin = 0
ymax = 0
zmin = 0.0
zmax = 2*d

nx = 100
nz = 100

dx = (xmax - xmin) / max(nx-1.0,1.0)
dz = (zmax - zmin) / max(nz-1.0,1.0)
xs = (np.r_[0:nx] - (nx-1.0)/2.0) * dx
zs = (np.r_[0:nz]) * dz

k = (2*np.pi)/lamda

factor = int(height / width)
ndiv = 32
#factor = 1

xs,zs = np.meshgrid(xs,zs,indexing='ij')

a = linear_array(nElements=nelex,pitch=width+kerf,kerf=kerf,height=height,c=soundspeed, nAbcissa=[ndiv,ndiv*factor])
a.focus([0,0,d],f0)

hmm = a.cw_pressure(xs,np.zeros(xs.shape),zs,k)
result = np.real(np.abs(hmm))
plt.figure()
plt.imshow(result,extent=np.round(100*np.r_[0,2*d,-d/2,d/2]),interpolation='none')
plt.xlabel('Depth [cm]')
plt.ylabel('Width [cm]')

