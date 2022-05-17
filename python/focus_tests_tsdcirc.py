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

plt.ion()

if __name__ == "__main__":

  eng = matlab.engine.start_matlab()
  focusPath = '/home/jmh/programming/matlab/focus'
  if not(os.path.exists(focusPath)):
    focusPath = 'c:/Users/jem/Documents/MATLAB/Focus'
  eng.addpath(eng.genpath('%s/Core' % (focusPath)))
  eng.chdir('%s/Examples' % (focusPath))

  print("Running Focus example...")
  eng.eval("run tsdcirc.m;",nargout=0)
  eng.eval("pressure = squeeze(pressure);",nargout=0)
  eng.eval("soundspeed = lossless.soundspeed;",nargout=0)

  print("Extracting input and output...")
  x = eng.workspace["x"]
  z = eng.workspace["z"]
  fs = eng.workspace["fs"]
  radius  = eng.workspace["radius"]
  f0      = eng.workspace["f0"]
  ncycles = eng.workspace["ncycles"]
  pressure = eng.workspace["pressure"]
  soundspeed = eng.workspace["soundspeed"]
  ndiv = int(eng.workspace["ndiv"])
  nx = int(eng.workspace["nx"])
  nz = int(eng.workspace["nz"])

  # Focus example uses (nx,nz) but gets (nx+1,nz+1) points
  nx = nx + 1
  nz = nz + 1

  d = radius * 2

  # Faster to supply info about Fortran ordering
if os.name == 'nt':
  reference = np.array(pressure).reshape(pressure.size, order='F')
else:
  reference = np.array(pressure._data).reshape(pressure.size, order='F')

xs,zs = np.meshgrid(x,z,indexing='ij')

print("Starting Python Engine...")
import addpaths
import swig_fnm as fnm

a = fnm.CircularApertureFloat(radius)
a.f0 = f0
a.nDivA = ndiv
a.fs = fs
a.w = ncycles / a.f0

pos = np.c_[xs.flatten(), 0.0*np.ones(nx*nz), zs.flatten()].astype(np.float32)

print("Running Sofus example...")
start = timer()
#bum = a.CalcTransientFieldRef(pos)[1].astype(np.float64) # Works
bum = a.CalcPwFieldThreaded(pos)[1].astype(np.float64)
end = timer()
tSpent = end - start
print('\nSeconds elapsed: %f' % (tSpent))
pressure = bum.reshape((nx, nz, bum.shape[1]))
pressure = pressure[:,:,0:reference.shape[2]]

# TODO: Consider tighting up BoxTimes - 3 samples are added (rounding, range, far field)
fh = plt.figure()
axes = [fh.add_subplot(121), fh.add_subplot(122)]

vpmin = pressure.min()
vpmax = pressure.max()

nFrames = pressure.shape[2]
pdata = []
for i in range(nFrames):
  pdata.append(pressure[:,:,i])

axes[0].imshow(reference[:,:,25])
axes[0].set_title('Focus reference')
imp = axes[1].imshow(pdata[0], vmin=vpmin, vmax=vpmax, animated=True)
axes[1].set_title('SOFUS result')

j = 0
def updatepfig(*args):
  global j, pdata;
  j = (j + 1) % nFrames
  imp.set_array(pdata[j])
  return imp,

pani = animation.FuncAnimation(fh, updatepfig,
                               interval=20, blit=True)

pressure = pressure[:,:,0:reference.shape[2]]
diff = np.fabs(pressure - reference)

denom = 0.5 * (np.fabs(pressure) + np.fabs(reference))
denom = np.maximum(denom,np.finfo(np.float64).eps)
# diff < -80 dB ignored
#diff[np.where(diff < 0.0001 * denom.max())] = 0.0
diff1 = 100.0 * diff / denom

print("Relative errors: %f (95%% percentile)" % (np.percentile(diff1.flatten(),95)))

data = []
nFrames = diff1.shape[2]
for i in range(nFrames):
  data.append(diff1[:,:,i])

cmap = plt.get_cmap('gist_rainbow')

vmin = diff1.min()
vmax = diff1.max()

fig = plt.figure()
im = plt.imshow(data[0], aspect='auto', cmap=cmap, extent=np.round(1000*np.r_[0,d,1.5*d/2,-1.5*d/2]), vmin=vmin, vmax=vmax, animated=True)
plt.colorbar()

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
ani = animation.FuncAnimation(fig, updatefig,
                              interval=20, blit=True)
eng.quit()

# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
