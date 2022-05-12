# Issue is when nDivH is odd and no offset in y (v)
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

swap = False#True

show3D = True
show3D = False

#cmap = plt.get_cmap('gist_rainbow')
cmap = parula_map

print('=================================[ transientrect.m ]=================================\n\n');
print('This example calculates the transient pressure field of a rectangular piston using a\n');
print('transient version of the Fast Nearfield Method. It outputs an animation of the\n');
print('pressure through time.\n\n');

width = 15e-3
height = 15e-3
d = width

f0=1e6
soundspeed = 1500 # m/s
lamda = soundspeed / f0 # wavelength, m

scale = 1.5

# Coordinate grid parameters
xmin = -scale*d/2
xmax = scale*d/2
ymin = 0
ymax = 0
zmin = 0
zmax = d

nx = 100
nz = 100
nx = 130
nz = 250

# To match that of Focus
if 1:
  nx = nx + 1
  nz = nz + 1

dx = (xmax - xmin)/nx
dz = (zmax - zmin)/nz

fs = 10e6
deltat = 1/fs

xs = (np.r_[0:nx] - (nx-1.0)/2.0) * dx
zs = (np.r_[0:nz]) * dz

k = (2*np.pi)/lamda

xs,zs = np.meshgrid(xs,zs,indexing='ij')

if 1:
  ys = 0.0*np.ones(xs.shape)
else:
  ys = -lamda/2.0*np.ones(xs.shape)

pos = np.c_[xs.flatten(), ys.flatten(), zs.flatten()].astype(np.float32)
ndiv = 40

a = fnm.ApertureFloat(1, width, 0.0, height)

a.f0 = f0
a.c = soundspeed
a.nDivH = ndiv
a.nDivW = ndiv
a.focus_type = 1

# Calculate the pressure
a.fs = fs

nCycles = 2*1.5

a.w = nCycles / a.f0

from timeit import default_timer as timer

ztest =  8e-3
xtest = 7.14e-3

start = timer()

a.nthreads = 4

pbar = fnm.ProgressBarStdOut()
a.ProgressBarSet(pbar)

# Not working - pulse gets to short (for positive delays)
#               no impact for negative delays
#a.delays = [30.0/fs]
bum = a.CalcTransientSingleElementNoDelay(pos,0x1F)[1]

end = timer()
tSpent = end - start
print('\nSeconds elapsed: %f' % (tSpent))

tSpent64 = tSpent * 10000.0 / (nx*nz) * 64
print('\nSeconds for 100x100 pixels, 64 channels: %f' % (tSpent64))

tSpentOptimized = (tSpent64 / 1.0)
print('\nSeconds for 100x100 pixels, 64 channels: %f (threaded)' % (tSpentOptimized))



# animation test

from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation

bum0 = bum.T
bum = bum0.reshape((bum0.shape[0], nx, nz))

vmin = bum.min()
vmax = bum.max()

data = []
nFrames = bum.shape[0]
for i in range(bum.shape[0]):
  data.append(bum[i,:,:])

if show3D:
  fig = plt.figure()
  ax = axes3d.Axes3D(fig)
  from matplotlib import cm

  wframe = ax.plot_wireframe(xs, zs, data[0], rstride=5, cstride=5)
  ax.set_zlim(vmin, vmax)

  # colors = cm.viridis(data[0])
  # rcount, ccount, _ = colors.shape
  # wframe = ax.plot_surface(xs, zs, data[0], rcount=rcount, ccount=ccount,
  #                      facecolors=colors, shade=False)
  # wframe.set_facecolor((0,0,0,0))
  # fig.show()

  def update(i, ax, fig, data):
      ax.cla()
      wframe = ax.plot_wireframe(xs, zs, data[i], rstride=5, cstride=5)
      ax.set_zlim(vmin,vmax)
      #colors = cm.viridis(data[i])
      #rcount, ccount, _ = colors.shape
      #wframe = ax.plot_surface(xs, zs, data[i], rcount=rcount, ccount=ccount,
      #                         facecolors=colors, shade=False)
      #wframe.set_facecolor((0,0,0,0))
      #fig.show()

      return wframe,

  ani = animation.FuncAnimation(fig, update,
                                frames=xrange(nFrames),
                                fargs=(ax, fig, data), interval=2)
else:
  fig = plt.figure()
  im = plt.imshow(data[0], aspect='auto', cmap=cmap, extent=np.round(1000*np.r_[0,d,1.5*d/2,-1.5*d/2]), vmin=vmin, vmax=vmax, animated=True)

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
