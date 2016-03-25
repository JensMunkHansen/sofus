import sys
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import addpaths
import swig_fnm as fnm
#from figtools import *
#from parula import parula_map

plt.ion()

def log_compress(pressure,dBrange=60):
  logp = np.abs(pressure)
  logp = logp / logp.max()
  logp = 20*np.log10(logp)
  if (dBrange != None):
      logp[logp < -dBrange] = -dBrange
      logp[0,0] = 0
      logp[-1,-1] = -dBrange
  return logp

f0 = 1e6
soundspeed = 1500
lambd = soundspeed / f0
width = 5 * lambd
height = 7.5 * lambd
nElements = 1
fs = 100e6

hh = height / 2.0
hw = width / 2.0
a = fnm.ApertureFloat(1,width, 0.0, height)

a.fs = fs
a.c  = soundspeed
a.f0 = f0

nx = 61
nz = 101 # Reduce from 81x101 to 81x86

sym = False
# Why do we need symmetry of x (need for compact and ultra)

if sym:
    xs = np.linspace(-1.5 * width/2, 1.5 * width/2, nx)
else:
    xs = np.linspace(0, 1.5 * width/2, nx)

zs = np.linspace(0,(width/2.0)**2 / lambd, nz) # 1 cm

[xs,zs] = np.meshgrid(xs,zs)
xs = xs.flatten()
zs = zs.flatten()
pos = np.c_[xs,np.zeros(nx*nz),zs].astype(np.float32)
xs = xs.reshape((nz,nx))
zs = zs.reshape((nz,nx))

a.method = 0x14
a.delays = np.zeros(1,dtype=np.float32)

out = a.CalcCwFast(pos)

out = out.reshape((nz,nx))

rho = 10000
omega = 2*np.pi*f0
pressure = -1j * omega * rho * np.exp(1j * omega * 0) * out

logp = log_compress(abs(pressure),None)
plt.figure()
plt.imshow(logp)
sys.exit(0)

xs = 1.5 * (xs / xs.max())
zs = 1.0 * (zs / zs.max())

latexify(columns=1)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xs,zs,logp, cmap=parula_map,rstride=2, cstride=2,edgecolors='k')
ax.view_init(elev=2*22.5, azim=-135)
ax.invert_xaxis()
ax.set_xlabel('Lateral distance [a]')
ax.set_ylabel('Axial distance [$a^2/\lambda$]')
ax.set_zlabel('Normalized pressure [dB]',rotation=180)
plt.tight_layout()
plt.show()
plt.savefig("./p_lockwood_rect.pdf", bbox_inches='tight')
