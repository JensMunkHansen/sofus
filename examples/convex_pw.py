import matplotlib.pyplot as plt
import numpy as np

import addpaths
import swig_fnm as fnm
from grid import grid

plt.ion()

nElements = 160
width = 0.32e-3
kerf = 1.5e-5
pitch = 0.33e-3
radius = 45.28e-3
height = 1.35e-2
efocus = 7e-2
density = 919.60 # kg/m^3
c = 1540.0
fxmt = 6.0
f0 = 2e6
fs = 5*f0

nDivH = int(18)
nDivW = int(4)

a = fnm.ApertureFloat_FocusedConvexArrayCreate(nElements,
                                               width,
                                               kerf,
                                               height,
                                               radius,
                                               1,
                                               0.0)[1]
a.f0 = f0
a.c  = c
a.fc = f0
a.fs = fs
a.rho = density
a.nDivH = nDivH
a.nDivW = nDivW

# Attenuation
#a.alpha =
#a.beta  =
a.att_enabled = False

# Set focusing type to Pythagorean
a.focus_type = fnm.FocusingType.Pythagorean
# This will update the delays (when needed)
a.focus = [0.0, 0.0, 9.65e-2]

# Set apodization type to rectangular (will force computation of apodization)
a.apodization_type = fnm.ApodizationType.ApodizationTypeRectangular
# Apodization set using f-number
a.XmtFNumberSet(fxmt)

# Length of pulse
a.w = 2 * (1.0/f0)

a.excitation_type = fnm.ExcitationType.ExcitationTypeToneBurst

# Length of pulse
a.w = 2 * (1.0/f0)


# Not very elegant
myGrid = grid(nx=100, dx=0.4e-3, nz=100, offset_z=51, dz=1.38e-3)
pos = myGrid.values()

# Compute pulsed-wave field
myField = a.CalcPwFnmThreaded(pos)[1]
img = np.sum(myField**2 / fs,axis=1).reshape((100,100))

logimg = 20*np.log10(img/img.max())

plt.figure()
plt.imshow(logimg,extent=[pos[:,2].min(), pos[:,2].max(), pos[:,0].min(), pos[:,0].max()])
