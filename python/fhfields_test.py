import sys, os
import numpy as np

import matplotlib.pyplot as plt

filedir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(filedir, '../fnm'))

from dicts import dotdict
from fhfields import (FraunhoferAperture, FraunhoferField)
from parula import parula_map

plt.ion()

opt = dotdict({'pitch0' : 0.2e-3,
               'N0'     : 64,
               'pitch1' : 0.2e-3,
               'N1'     : 64,
               'f0'     : 7e6,
               'focus'  : np.r_[0.0,0.0,3.0e-2]})

xdc = FraunhoferAperture(**opt)
xdc.bandwidth = 0.2

resolution = np.r_[1.000e-004, 1.000e-004, 1.0000e-005]
rng        = np.r_[0.003, 0.003, 0.006]

opt = dotdict({'resolution' : resolution,
               'rng'        : rng})

[P12, fx1, fy1, f1] = FraunhoferField.FieldAtFocusNorway(xdc, **opt)

[p12,xax,yax,zax] = FraunhoferField.FourierToSpaceTimeNorway(P12,fx1,fy1,f1,**opt)

# p12 is (nx,ny,nz)
p1xy = 10*np.log10(np.mean(np.abs(p12)**2.0,axis=2)).T # (ny,nx)
p1xz = 20*np.log10(np.maximum(0, np.abs(p12[:,int(p12.shape[1]/2),:]))).T # (nz,nx)

x=1e3*xax
y=1e3*yax
z=1e3*zax;

dBrange = 50

fh = plt.figure()
axes = [fh.add_subplot(321),
        fh.add_subplot(323),
        fh.add_subplot(325),
        fh.add_subplot(322),
        fh.add_subplot(324),
        fh.add_subplot(326)]
p1xys = p1xy - p1xy.max()
im = axes[0].imshow(p1xys, extent=[x.min(), x.max(), y.min(), y.max()], cmap=parula_map)
im.set_clim(-dBrange,0)
p1xzs = p1xz - p1xz.max()
im = axes[1].imshow(p1xzs, extent=[x.min(), x.max(), z.min(), z.max()], cmap=parula_map)
im.set_clim(-dBrange,0)
axes[2].plot(p1xy[int(p1xy.shape[0]/2),:] - p1xy.max())
axes[2].set_ylim([-dBrange, 0])
