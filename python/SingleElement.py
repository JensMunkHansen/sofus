import numpy as np

from fnm import rect

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

nx = 50
nz = 70

dx = (xmax - xmin) / max(nx-1.0,1.0)
dz = (zmax - zmin) / max(nz-1.0,1.0)
xs = (np.r_[0:nx] - (nx-1.0)/2.0) * dx
zs = (np.r_[0:nz]) * dz

ndiv = 4
rect = rect(hw=width/2.0,hh=height/2.0,nAbcissa=[4,4])

k = (2*np.pi)/lamda

start = timer()
xs,zs = np.meshgrid(xs,zs,indexing='ij')

result2 = rect.H4(xs,np.zeros(xs.shape),zs,k)

end = timer()
timeString = create_time_string(end-start)
print(timeString)

plt.figure()
result2 = np.real(np.abs(result2))
plt.imshow(result2,extent=np.round(1000*np.r_[0,2*d,-d/2,d/2]),interpolation='none')
plt.xlabel('Depth [mm]')
plt.ylabel('Width [mm]')
