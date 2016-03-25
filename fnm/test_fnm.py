import addpaths
import numpy as np
import sys
import swig_fnm as fnm

import numpy as np

from fnm import rect

from timeit import default_timer as timer

import matplotlib.pyplot as plt
plt.ion()

def log_compress(pressure,dBrange=20):
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

f0 = 1e6 # excitation frequency,Hz
soundspeed = 1500 # m/s
lamda = soundspeed / f0 # wavelength, m

#define a transducer structure/array
nelex = 1
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

# Fast if many z-coordinates
nx = 130
nz = 250

dx = (xmax - xmin) / max(nx-1.0,1.0)
dz = (zmax - zmin) / max(nz-1.0,1.0)
xs = (np.r_[0:nx] - (nx-1.0)/2.0) * dx
zs = (np.r_[0:nz]) * dz

factor = int(height / width) # 16
ndiv = 4

k = (2*np.pi)/lamda

xs,zs = np.meshgrid(xs,zs,indexing='ij')

ys = np.zeros(xs.shape)

rect = rect(hw=width/2.0,hh=height/2.0,nAbcissa=[ndiv,ndiv*factor])

if 1:
  start = timer()
  result2 = rect.H4(xs,ys,zs,k)
  end = timer()

  timeString = create_time_string(end-start)
  print(timeString)

  plt.figure()
  result2 = np.real(np.abs(result2))
  plt.imshow(log_compress(result2),aspect='auto',extent=np.round(1000*np.r_[0,2*d,-d/2,d/2]),interpolation='none')
  plt.xlabel('Depth [mm]')
  plt.ylabel('Width [mm]')

if 1:
  start = timer()
  result2 = rect.HN(xs,ys,zs,k)
  end = timer()

  timeString = create_time_string(end-start)
  print(timeString)

  plt.figure()
  result2 = np.real(np.abs(result2))
  result2 = result2.reshape((nx,nz))
  plt.imshow(log_compress(result2),aspect='auto',extent=np.round(1000*np.r_[0,2*d,-d/2,d/2]),interpolation='none')
  plt.xlabel('Depth [mm]')
  plt.ylabel('Width [mm]')
  

pos = np.c_[xs.flatten(), ys.flatten(), zs.flatten()].astype(np.float32)

a = fnm.ApertureFloat(1,width,kerf,height)
a.f0 = f0
a.c  = soundspeed
a.nDivH = ndiv
a.nDivW = ndiv*factor

start = timer()
out = a.CalcCwFast(pos)
end = timer()
timeString = create_time_string(end-start)
print(timeString)

plt.figure()
result3 = np.real(np.abs(out)).reshape((nx,nz))
plt.imshow(log_compress(result3),aspect='auto',extent=np.round(1000*np.r_[0,2*d,-d/2,d/2]),interpolation='none')
plt.xlabel('Depth [mm]')
plt.ylabel('Width [mm]')


