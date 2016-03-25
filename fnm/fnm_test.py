import sys
import numpy as np
import fnm as ref

import addpaths
import swig_fnm as lib

import matplotlib.pyplot as plt
plt.ion()

ndiv = 6

r = ref.rect(hh=0.5,hw=0.5,nAbcissa=[ndiv,ndiv])

point = [0,0,5.0]
k = 1.0

dirt = r.H_ref(point,k)

fast = r.H(np.ones((1,1))*point[0],
           np.ones((1,1))*point[1],
           np.ones((1,1))*point[2],
           k)

assert(np.abs(dirt-fast[0,0]) < 1e-6)

acc = r.H4(np.ones((1,1))*point[0],
            np.ones((1,1))*point[1],
            np.ones((1,1))*point[2],
            k)


a = lib.ApertureFloat(1,1.0,0.0,1.0)
a.nDivW = ndiv
a.nDivH = ndiv
a.c = 2*np.pi
a.f0 = 1.0
a.nthreads = 4
a.focus = [0,0,d]

pos = np.r_[[point]].astype(np.float32)
out  = a.CalcCwFieldRef(pos)

nx = 200
dx = 0.05
xs = (np.r_[0:nx] - (nx-1.0)/2.0) * dx
ys = 0.0*np.ones(nx)
pos  = np.c_[xs, ys, point[2]*np.ones(nx)].astype(np.float32)
out1 = a.CalcCwFieldRef(pos)
out2 = a.CalcCwFast(pos)

out1acc = r.H4(np.c_[pos[:,0]],np.c_[pos[:,1]],np.c_[pos[:,2]],k)

# Looks better
out1fast = r.H(np.c_[pos[:,0]],np.c_[pos[:,1]],np.c_[pos[:,2]],k)

out2fast = r.HN(np.c_[pos[:,0]],np.c_[pos[:,1]],np.c_[pos[:,2]],k)

plt.figure()
plt.plot(np.abs(out2),'k')

out1ref = np.r_[[r.H_ref(pos[i,:],k) for i in range(nx)]] # Divided with 2*pi
plt.plot(np.abs(out1acc),'b')
plt.plot(np.abs(out1ref),'r')
plt.plot(np.abs(out1),'g')
plt.plot(np.abs(out1fast),'m')
plt.plot(np.abs(out2fast),'y')
plt.legend(['SIMD', 'H4', 'H_ref', 'C++', 'PyFast', 'HN'])
