# TODO: Make double precision work

import os
import sys
import numpy as np
try:
  import matlab.engine
except ModuleNotFoundError:
  sys.exit(-1)
import matplotlib.pyplot as plt

from timeit import default_timer as timer

import matplotlib.animation as animation

doublePrecision = False

global i, j, k
i = 0
j = 0
k = 0

if (doublePrecision):
  dtype = np.float64
else:
  dtype = np.float32

print("Starting Matlab Engine...")
eng = matlab.engine.start_matlab()
focusPath = '/home/jmh/programming/matlab/focus'
if not(os.path.exists(focusPath)):
  focusPath = 'c:/Users/jem/Documents/MATLAB/Focus'

eng.addpath(eng.genpath('%s/Core' % (focusPath)))
eng.chdir('%s/Examples' % (focusPath))

print("Running Focus example...")
eng.eval("run tsdrectarray.m;",nargout=0)
eng.eval("pressure = squeeze(pressure);",nargout=0)
eng.eval("soundspeed = lossless.soundspeed;",nargout=0)

print("Extracting input and output...")
x = eng.workspace["x"]
z = eng.workspace["z"]
fs = eng.workspace["fs"]
width  = eng.workspace["width"]
kerf = eng.workspace["kerf"]
height = eng.workspace["height"]
nelex = int(eng.workspace["nelex"])
neley = int(eng.workspace["neley"])

f0      = eng.workspace["f0"]
ncycles = eng.workspace["ncycles"]
pressure = eng.workspace["pressure"]
soundspeed = eng.workspace["soundspeed"]
ndiv = int(eng.workspace["ndiv"])
nx = int(eng.workspace["nx"])
nz = int(eng.workspace["nz"])

d = nelex * (width+kerf);

# Focus example uses (nx,nz) but gets (nx+1,nz+1) points
nx = nx + 1
nz = nz + 1

ndiv = ndiv * 4

# Faster to supply info about Fortran ordering
if os.name == 'nt':
  reference = np.array(pressure).reshape(pressure.size, order='F')
else:
  reference = np.array(pressure._data).reshape(pressure.size, order='F')

# Testing remove near field
#z = np.array(z)
#z = z[0,50:]
#nz = len(z)

xs,zs = np.meshgrid(x,z,indexing='ij')

print("Starting Python Engine...")
import addpaths
import swig_fnm as fnm

if (doublePrecision):
  a = fnm.ApertureDouble(nelex, width, kerf, height)
else:
  a = fnm.ApertureFloat(nelex, width, kerf, height)

a.f0 = f0
a.nDivH = ndiv
a.nDivW = ndiv
a.focus_type = fnm.FocusingType.Pythagorean
a.fs = fs
a.c = soundspeed
a.w = ncycles / a.f0
a.focus = [0,0,d]

a.FocusUpdate()
delays = a.delays
a.focus_type = fnm.FocusingType.Delays
delays = delays - delays.min()
a.delays = delays

a.nthreads = 4

a.excitation_type = fnm.ExcitationType.ExcitationTypeToneBurst

pbar = fnm.ProgressBarStdOut()
a.ProgressBarSet(pbar)

pos = np.c_[xs.flatten(), 0.0*np.ones(nx*nz), zs.flatten()].astype(dtype)

ix = 75
iz = 100
iPoint = ix*nz + iz
#pos = pos[[iPoint],:]

print("Running Sofus example...")
start = timer()
bum = a.CalcPwFnmThreaded(pos,0x1F)[1]
end = timer()
tSpent = end - start
print('\nSeconds elapsed: %f' % (tSpent))
pressure = bum.reshape((nx, nz, bum.shape[1]))
pressure = pressure[:,:,0:reference.shape[2]]

# TODO: Consider tighting up BoxTimes - 3 samples are added (rounding, range, far field)
fh0 = plt.figure()
axes = [fh0.add_subplot(121), fh0.add_subplot(122)]

vpmin = pressure.min()
vpmax = pressure.max()

nFrames = pressure.shape[2]
pdata = []
for i in range(nFrames):
  pdata.append(pressure[:,:,i])

vmin = reference.min()
vmax = reference.max()

data = []
for i in range(nFrames):
  data.append(reference[:,:,i])

im = axes[0].imshow(data[0], vmin=vmin, vmax=vmax, animated=True)
axes[0].set_title('Focus reference')
imp = axes[1].imshow(pdata[0], vmin=vpmin, vmax=vpmax, animated=True)
axes[1].set_title('SOFUS result')

j = 0
def updatepfig(*args):
  global j, pdata;
  j = (j + 1) % nFrames
  imp.set_array(pdata[j])
  return imp,

pani = animation.FuncAnimation(fh0, updatepfig,
                               interval=20, blit=True)

i = 0
def updatefig(*args):
  global i, data;
  i = (i + 1) % nFrames
  im.set_array(data[i])
  return im,

ani = animation.FuncAnimation(fh0, updatefig,
                              interval=20, blit=True)

plt.show()

fh1 = plt.figure()
axes = [fh1.add_subplot(121), fh1.add_subplot(122)]
refimg = np.sum(reference**2,axis=2)
logref = 20*np.log10(refimg/refimg.max())
img = axes[0].imshow(logref)
img.set_clim([-60, 0])
axes[0].set_title('Focus reference')
pimg = np.sum(pressure**2,axis=2)
logp = 20*np.log10(pimg/pimg.max())
img = axes[1].imshow(logp)
img.set_clim([-60, 0])
axes[1].set_title('SOFUS result')
plt.show()

diff = np.fabs(pressure - reference)

denom = 0.5 * (np.fabs(pressure) + np.fabs(reference))
denom = np.maximum(denom,np.finfo(np.float32).eps)
# diff < -60 dB ignored
diff[np.where(diff < 0.001 * denom.max())] = 0.0
diff1 = 100.0 * diff / denom

print("Relative errors: %f (95%% percentile)" % (np.percentile(diff1.flatten(),95)))

ddata = []
nFrames = diff1.shape[2]
for i in range(nFrames):
  ddata.append(diff1[:,:,i])

cmap = plt.get_cmap('gist_rainbow')

vdmin = diff1.min()
vdmax = diff1.max()

fig = plt.figure()
imd = plt.imshow(ddata[0], aspect='auto', cmap=cmap, extent=np.round(1000*np.r_[0,d,1.5*d/2,-1.5*d/2]), vmin=vdmin, vmax=vdmax, animated=True)
plt.colorbar()

k = 0
def updatedfig(*args):
  global k, data;
  k = (k + 1) % nFrames
  imd.set_array(ddata[k])
  return imd,

anim_running = True

def onClick(event):
  global anim_running
  if anim_running:
    anid.event_source.stop()
    anim_running = False
  else:
    anid.event_source.start()
    anim_running = True

fig.canvas.mpl_connect('button_press_event', onClick)
anid = animation.FuncAnimation(fig, updatedfig,
                               interval=20, blit=True)

np.savez('issue.npz', reference=reference[ix,iz,:], pressure=pressure[ix,iz,:])
import pickle
pickle.dump(a,open('xdc.p','wb'))

plt.show()
eng.quit()

# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
