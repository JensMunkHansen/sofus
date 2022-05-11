import sys
import numpy as np

import matplotlib.pyplot as plt

import addpaths
import swig_fnm as fnm

from parula import parula_map

def log_compress(pressure,dBrange=60):
  logp = np.abs(pressure)
  logp = logp / logp.max()
  logp = 20*np.log10(logp)
  logp[logp < -dBrange] = -dBrange
  logp[0,0] = 0
  logp[-1,-1] = -dBrange
  return logp

plt.ion()

cmap = plt.get_cmap('gist_rainbow')
cmap = parula_map

print('============================[ test_circ_transient.py ]============================\n\n')
print('This example calculates the transient pressure field of a circular piston using\n')
print('the Fast Nearfield Method. \n\n')

# demo file for FNM transient.
# define a xdcr structure/array
radius = 7.5e-3
d = radius * 2

# center frequency 1MHz
f0=1e6
soundspeed = 1500 # m/s
lamda = soundspeed / f0 # wavelength, m

# define the calculation grid
xmin = -1.5*d/2
xmax = 1.5*d/2
ymin = 0
ymax = 0
zmin = 0
zmax = d

nx = 150
nz = 300

dx = (xmax - xmin) / max(nx-1.0,1.0)
dz = (zmax - zmin) / max(nz-1.0,1.0)
#xs = (np.r_[0:nx] - (nx-1.0)/2.0 +0.5) * dx # TODO: Fix x = 0.0
xs = (np.r_[0:nx] - (nx-1.0)/2.0) * dx
zs = (np.r_[0:nz] + 0.5) * dz

ndiv=20

k = (2*np.pi)/lamda

xs,zs = np.meshgrid(xs,zs,indexing='ij')

ys = 0.0*np.ones(xs.shape)

pos = np.c_[xs.flatten(), np.zeros(nx*nz), zs.flatten()].astype(np.float32)

a = fnm.CircularApertureFloat(radius)
a.f0 = f0
a.nDivA = ndiv
a.gridSectorScale = 1.0

# run the function!
pressure = a.CalcCwFieldRef(pos)[1]


pressure = pressure.reshape((nx,nz))

pressure = np.abs(pressure)
pressure = pressure.astype(np.float64)
logp = log_compress(pressure,30)

fig = plt.figure()
ax = fig.add_subplot(121)

ax.imshow(logp,aspect='auto',
          cmap=cmap,
          extent=np.round(1000*np.r_[0,d,-1.5*d/2,1.5*d/2]),interpolation='none')

plt.title('Pressure at y = 0')
plt.xlabel('z (mm)')
plt.ylabel('x (mm)')

a.fs = 10e6

nCycles = 3.0

a.w = nCycles / a.f0

from timeit import default_timer as timer

#bum = a.CalcTransientFieldRef([pos[0,:]])[1]

#sys.exit(0)

start = timer()
bum = a.CalcTransientFieldRef(pos)[1]
end = timer()
print('\nSeconds elapsed: %f' % (end - start))

ax = fig.add_subplot(122)
p = np.sum(np.abs(bum**2),axis=1).reshape((nx,nz))
ax.imshow(p,aspect='auto',
          cmap=cmap,
          extent=np.round(1000*np.r_[0,d,-1.5*d/2,1.5*d/2]),interpolation='none')

# animation test

from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation


bum = bum.T
bum = bum.reshape((bum.shape[0], nx, nz))

data = []
nFrames = bum.shape[0]
for i in range(bum.shape[0]):
  data.append(bum[i,:,:])

vmin = bum.min()
vmax = bum.max()
  
if 0:
  fig = plt.figure()
  ax = axes3d.Axes3D(fig)
  
  wframe = ax.plot_wireframe(xs, zs, data[0], rstride=2, cstride=2)
  ax.set_zlim(vmin, vmax)
  
  def update(i, ax, fig, data):
      ax.cla()
      wframe = ax.plot_wireframe(xs, zs, data[i], rstride=5, cstride=5)
      ax.set_zlim(vmin,vmax)
      return wframe,
  
  ani = animation.FuncAnimation(fig, update, 
                                frames=xrange(nFrames),
                                fargs=(ax, fig, data), interval=2)
  plt.show()

if 1:
  fig = plt.figure()
  im = plt.imshow(data[0], aspect='auto', cmap=cmap, extent=np.round(1000*np.r_[0,d,-1.5*d/2,1.5*d/2]), vmin=vmin, vmax=vmax, animated=True)

  i = 0
  def updatefig(*args):
    global i, data;
    i = (i + 1) % nFrames
    im.set_array(data[i])
    return im,

  anim_running = True

  def onClick(event):
    global anim_running
    if anim_running:
      ani.event_source.stop()
      anim_running = False
    else:
      ani.event_source.start()
      anim_running = True  
  
  fig.canvas.mpl_connect('button_press_event', onClick)
  # Could we add two functions, simPoints and simData
  ani = animation.FuncAnimation(fig, updatefig, 
                                interval=20, blit=True)
  plt.show()

