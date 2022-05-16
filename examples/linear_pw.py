import matplotlib.pyplot as plt
import numpy as np

import addpaths
import swig_fnm as fnm
from grid import grid

plt.ion()

nElements = 64
width = 0.18e-3
kerf = 2e-5
pitch = 2.0e-4
height = 1.2e-2
density = 919.60 # kg/m^3
c = 1540.0
fxmt = 2.0
f0 = 7e6
fs = 5*f0

nDivH = int(56)
nDivW = int(4)

a = fnm.ApertureFloat_FocusedLinearArrayCreate(nElements,
                                               width,
                                               kerf,
                                               height,
                                               1,
                                               0.0)[1]

a.c  = c
a.fs = fs

nCycles = 2

#a.excitation_type = fnm.ExcitationType.ExcitationTypeToneBurst # Error here
a.excitation_type = fnm.ExcitationType.ExcitationTypeHanningWeightedPulse

# Length of pulse
a.w = nCycles * (1.0/f0)

a.f0 = f0
a.nDivH = nDivH
a.nDivW = nDivW
a.focus = [0.0, 0.0, 3e-2]
a.rho = density

# Set focusing type to Pythagorean
a.focus_type = fnm.FocusingType.Pythagorean

a.fc = f0

# Attenuation
#a.alpha =
#a.beta  =
a.att_enabled = False

# This will update the delays (when needed)

# Set apodization type to rectangular (will force computation of apodization)
a.apodization_type = fnm.ApodizationType.ApodizationTypeRectangular
# Apodization set using f-number
a.XmtFNumberSet(fxmt)

a.nthreads = 8

# Not very elegant
myGrid = grid(nx=100, dx=6e-5, nz=100, offset_z=60, dz=4e-4)
pos = myGrid.values()

# Compute pulsed-wave field
myField = a.CalcPwFnmThreaded(pos)[1]
img = np.sum(myField**2/fs,axis=1).reshape((100,100))

logimg = 20*np.log10(img/img.max())

plt.imshow(logimg)
