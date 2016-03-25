import addpaths
import numpy as np
import sys
import swig_fnm as fnm

import numpy as np
from scipy.signal import hilbert

from fnm import (rect,linear_array)

from timeit import default_timer as timer

import matplotlib.pyplot as plt
plt.ion()

def log_compress(pressure,dBrange=60):
  logp = np.abs(pressure)
  logp = logp / logp.max()
  logp = 20*np.log10(logp)
  logp[logp < -dBrange] = -dBrange
  logp[0,0] = 0
  logp[-1,-1] = -dBrange
  return logp
  
  

def create_time_string(seconds):
  m, s = divmod(seconds, 60)
  h, m = divmod(m, 60)
  ms = np.round(s * 1000) % 1000
  timeString = "%d:%02d:%02d,%03d" % (h, m, s, (1000*ms))
  return timeString;

print('This script uses the Fast Nearfield Method to calculate the CW pressure field of');
print('a single element. The script outputs the pressure field.\n');

start = timer()

f0 = 1e6 # excitation frequency,Hz
soundspeed = 1500 # m/s
lamda = soundspeed / f0 # wavelength, m

#define a transducer structure/array
nelex = 128
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

nx = 170
nz = 250

dx = (xmax - xmin) / max(nx-1.0,1.0)
dz = (zmax - zmin) / max(nz-1.0,1.0)
xs = (np.r_[0:nx] - (nx-1.0)/2.0) * dx + 0.5*dx 
zs = (np.r_[0:nz]) * dz

factor = int(height / width)
ndiv = 20
rect = rect(hw=width/2.0,hh=height/2.0,nAbcissa=[ndiv,ndiv])

k = (2*np.pi)/lamda

start = timer()
xs,zs = np.meshgrid(xs,zs,indexing='ij')

ys = 0.0*np.ones(xs.shape)

pos = np.c_[xs.flatten(), np.zeros(nx*nz), zs.flatten()].astype(np.float32)

# Original
a = fnm.ApertureFloat(nelex,width,kerf,height)
a.f0 = f0
a.c  = soundspeed
a.nDivW = ndiv
a.nDivH = ndiv
a.nthreads = 4

d = nelex * (width+kerf)
a.focus = [0,0,d]

out = a.CalcCwFast(pos)
#out = a.CalcCwFieldRef(pos)
end = timer()
timeString = create_time_string(end-start)
print(timeString)

plt.figure()
result3 = out.reshape((nx,nz))

rho = 10000
omega = 2*np.pi*f0
pressure = -1j * omega * rho * np.exp(1j * omega * 0) * result3
# -j omega rho v exp(jwt)

logp = np.abs(pressure)#log_compress(np.abs(pressure))

plt.imshow(logp,aspect='auto',extent=np.round(1000*np.r_[0,2*d,-1.5*d/2,1.5*d/2]),interpolation='none')
plt.xlabel('Depth [mm]')
plt.ylabel('Width [mm]')

if 0:
  delays = np.sum((a.positions - a.focus)**2,axis=1)/a.c
  a.phases = 2*np.pi*a.f0 * (delays - delays.max())
  out = a.CalcCwFieldRef(pos)
  result3 = out.reshape((nx,nz))
  pressure = -1j * omega * rho * np.exp(1j * omega * 0) * result3

  logp = log_compress(pressure)
  plt.figure()
  plt.imshow(logp,extent=np.round(1000*np.r_[0,2*d,-1.5*d/2,1.5*d/2]),interpolation='none')
  plt.xlabel('Depth [mm]')
  plt.ylabel('Width [mm]')

if 0:
  a1 = linear_array(nElements=nelex,pitch=width+kerf,kerf=kerf,height=height,c=soundspeed, nAbcissa=[a.nDivW,a.nDivH])
  a1.set_focus([0,0,d],f0)
  result2 = a1.cw_pressure(xs,np.zeros(xs.shape),zs,k)
  logp1 = log_compress(result2)
  plt.figure()
  plt.imshow(logp1,aspect='auto',extent=np.round(1000*np.r_[0,2*d,-1.5*d/2,1.5*d/2]),interpolation='none')

