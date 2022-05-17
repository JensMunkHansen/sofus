import os
import numpy as np
import matlab.engine
import matplotlib.pyplot as plt

from timeit import default_timer as timer

import matplotlib.animation as animation

plt.ion()

doublePrecision = False

if doublePrecision:
  dtype = np.float64
else:
  dtype = np.float32

eng = matlab.engine.start_matlab()
focusPath = '/home/jmh/programming/matlab/focus'
if not(os.path.exists(focusPath)):
  focusPath = 'c:/Users/jem/Documents/MATLAB/Focus'
eng.addpath(eng.genpath('%s/Core' % (focusPath)))
eng.chdir('%s/Examples' % (focusPath))

print("Running Focus example...")
eng.eval("run SingleFocusRectArray.m;",nargout=0)
eng.eval("pressure = squeeze(pressure);",nargout=0)
eng.eval("soundspeed = lossless.soundspeed;",nargout=0)
eng.eval("rho = lossless.density;",nargout=0)

print("Extracting input and output...")
width  = eng.workspace["width"]
kerf = eng.workspace["kerf"]
height = eng.workspace["height"]
nelex = int(eng.workspace["nelex"])
neley = int(eng.workspace["neley"])
x = eng.workspace["x"]
z = eng.workspace["z"]
f0      = eng.workspace["f0"]
rho     = eng.workspace["rho"]
pressure = eng.workspace["pressure"]
soundspeed = eng.workspace["soundspeed"]
ndiv = int(eng.workspace["ndiv"])
nx = int(eng.workspace["nx"])
nz = int(eng.workspace["nz"])

# Focus example uses (nx,nz) but gets (nx+1,nz+1) points
nx = nx + 1
nz = nz + 1

# Faster to supply info about Fortran ordering
reference = np.array(pressure).reshape(pressure.size, order='F')

xs,zs = np.meshgrid(x,z,indexing='ij')

print("Starting Python Engine...")
import addpaths
import swig_fnm as fnm

a = fnm.ApertureFloat(nelex, width, kerf, height)

a.f0 = f0
a.nDivH = ndiv
a.nDivW = ndiv
a.c = soundspeed
a.nthreads = 4

# Only works for D-term
#a.delays = [0.5/ a.f0]

pbar = fnm.ProgressBarStdOut()
#a.ProgressBarSet(pbar)

pos = np.c_[xs.flatten(), 0.0*np.ones(nx*nz), zs.flatten()].astype(dtype)

print("Running Sofus example...")
start = timer()
# Small error signal around nearest edge from the other edge wave (x=-7.788 mm, ix=20, iz=7, 0x08)
# Outside aperture, so there should be a little signal. Is McGough computing a derivative requiring
# two abcissas
d = nelex * (width+kerf)
a.focus = [0,0,d]
a.focus_type = fnm.FocusingType.Rayleigh
a.rho = rho
pressure = a.CalcCwFast(pos)[1]

pressure = pressure.reshape((nx,nz))

# Taken care of inside library
# omega = 2*np.pi*f0
# pressure = -1j * omega * rho * np.exp(1j * omega * 0) * pressure

pressure = -1j * pressure
pressure = np.abs(pressure)
pressure = pressure.astype(np.float64)

end = timer()
tSpent = end - start
print('\nSeconds elapsed: %f' % (tSpent))

fh0 = plt.figure()
axes = [fh0.add_subplot(121), fh0.add_subplot(122)]
extent = [xs.min(), xs.max(), zs.min(), zs.max()]
im = axes[0].imshow(np.abs(reference), extent=extent)
axes[0].set_title('Focus reference')
imp = axes[1].imshow(np.abs(pressure), extent=extent)
axes[1].set_title('SOFUS result')

# We ignore face
diff = np.abs(pressure) - np.abs(reference)
denom = 0.5 * (np.abs(pressure) + np.abs(reference))
denom = np.maximum(denom,np.finfo(np.float32).eps)
# diff < -60 dB ignored
diff[np.where(diff < 0.001 * denom.max())] = 0.0
diff1 = 100.0 * diff / denom

print("Relative errors: %f (95%% percentile)" % (np.percentile(diff1.flatten(),95)))


# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
