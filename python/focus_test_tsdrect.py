# -*- coding: utf-8; tab-width: 2; python-indent: 2; indent-tabs-mode: nil -*-

import sys, os
import numpy as np
import matlab.engine
import matplotlib.pyplot as plt

from timeit import default_timer as timer

import matplotlib.animation as animation

doublePrecision = False
if doublePrecision:
  dtype = np.float64
else:
  dtype = np.float32

plt.ion()

print("Starting Matlab Engine...")
eng = matlab.engine.start_matlab()
focusPath = '/home/jmh/programming/matlab/focus'
if not(os.path.exists(focusPath)):
  focusPath = 'c:/Users/jem/Documents/MATLAB/Focus'

eng.addpath(eng.genpath('%s/Core' % (focusPath)))
eng.chdir('%s/Examples' % (focusPath))

print("Running Focus example...")
eng.eval("run tsdrect.m;",nargout=0)
eng.eval("pressure = squeeze(pressure);",nargout=0)
eng.eval("soundspeed = lossless.soundspeed;",nargout=0)

print("Extracting input and output...")
x          = eng.workspace["x"]
z          = eng.workspace["z"]
fs         = eng.workspace["fs"]
width      = eng.workspace["width"]
height     = eng.workspace["height"]
f0         = eng.workspace["f0"]
ncycles    = eng.workspace["ncycles"]
pressure   = eng.workspace["pressure"]
soundspeed = eng.workspace["soundspeed"]
ndiv       = int(eng.workspace["ndiv"])
nx         = int(eng.workspace["nx"])
nz         = int(eng.workspace["nz"])

global i, j

# Focus example uses (nx,nz) but gets (nx+1,nz+1) points
nx = nx + 1
nz = nz + 1

d = width

# Faster to supply info about Fortran ordering
if os.name == 'nt':
  reference = np.array(pressure).reshape(pressure.size, order='F')
else:
  reference = np.array(pressure._data).reshape(pressure.size, order='F')

xs,zs = np.meshgrid(x,z,indexing='ij')

print("Starting Python Engine...")
import addpaths
import swig_fnm as fnm

if doublePrecision:
  a = fnm.ApertureDouble(1, width, 0.0, height)
else:
  a = fnm.ApertureFloat(1, width, 0.0, height)

a.f0 = f0
a.nDivH = ndiv
a.nDivW = ndiv
a.focus_type = fnm.FocusingType.Pythagorean
a.fs = fs
a.c = soundspeed
a.w = ncycles / a.f0
a.excitation_type = fnm.ExcitationType.ExcitationTypeToneBurst

a.nthreads = 4

# Only works for D-term
#a.delays = [0.5/ a.f0]

pbar = fnm.ProgressBarStdOut()
a.ProgressBarSet(pbar)

pos = np.c_[xs.flatten(), 0.0*np.ones(nx*nz), zs.flatten()].astype(dtype)

if 0:
  ix = 20
  iz = 7
  iPoint = ix*nz + iz
  pos = pos[[iPoint],:]

print("Running Sofus example...")
start = timer()
# Small error signal around nearest edge from the other edge wave (x=-7.788 mm, ix=20, iz=7, 0x08)
# Outside aperture, so there should be a little signal. Is McGough computing a derivative requiring
# two abcissas

tstart, bum = a.CalcPwFnmThreaded(pos,0x1F)
end = timer()
tSpent = end - start
print('\nSeconds elapsed: %f' % (tSpent))
pressure = bum.reshape((nx, nz, bum.shape[1]))
pressure = pressure[:,:,0:reference.shape[2]].astype(np.float64)

# TODO: Consider tighting up BoxTimes - 3 samples are added (rounding, range, far field)
fh = plt.figure()
axes = [fh.add_subplot(121), fh.add_subplot(122)]

vpmin = pressure.min()
vpmax = pressure.max()

nFrames = pressure.shape[2]
pdata = []
for i in range(nFrames):
  pdata.append(pressure[:,:,i])

extent=np.round(1000*np.r_[0,d,1.5*d/2,-1.5*d/2])

axes[0].imshow(reference[:,:,25],extent=extent, vmin=vpmin, vmax=vpmax)
axes[0].set_title('Focus reference')
imp = axes[1].imshow(pdata[0], extent=extent, vmin=vpmin, vmax=vpmax, animated=True)
axes[1].set_title('SOFUS result')

j = 0
def updatepfig(*args):
  global j, pdata;
  j = (j + 1) % nFrames
  imp.set_array(pdata[j])
  return imp,

pani = animation.FuncAnimation(fh, updatepfig,
                               interval=20, blit=True)

diff = np.fabs(pressure - reference)

# 95% of points less than 0.5% error (global)
denom0 = 0.5 * (np.fabs(pressure).max() + np.fabs(reference).max())
diff0 = 100.0 * diff / denom0
err95 = np.percentile(diff0.flatten(),95)
print("95%% percentile: %f %% error (rel. to global max)" % (err95))

# More fair comparison
denom = 0.5 * (np.fabs(pressure) + np.fabs(reference))
denom = np.maximum(denom,np.finfo(np.float64).eps)
# diff < -60 dB ignored
diff[np.where(diff < 0.001 * denom.max())] = 0.0
diff1 = 100.0 * diff / denom

err95 = np.percentile(diff1.flatten(),95)
print("95%% percentile: %f %% error" % (err95))


data = []
nFrames = diff1.shape[2]
for i in range(nFrames):
  data.append(diff1[:,:,i])

cmap = plt.get_cmap('gist_rainbow')

vmin = diff1.min()
vmax = diff1.max()

fig = plt.figure()
im = plt.imshow(data[0], aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, animated=True)
plt.colorbar()
plt.title("Relative difference")
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
eng.quit()
