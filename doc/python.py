##
# @page python Python Sample
#
# TODO: Far field approximation has an error
import sys
import numpy as np
import matplotlib.pyplot as plt

import addpaths

import PyField2 as f2

import swig_fnm as fnm

plt.ion()

fs = 100e6
c  = 1500.0
f2.field_init(-1)
f2.set_field('use_lines', 0)
f2.set_field('use_rectangles', 1)
f2.set_sampling(fs)
f2.set_field('c', c)

width  = 2.5e-3
height = 5.0e-3
kerf   = 0.0e-3
nElements = 2
xdc = f2.xdc_linear_array(nElements, width, height, kerf, 1, 1, [0,0,0])

impulse_response = np.r_[1]
excitation       = np.r_[1]

f2.xdc_impulse(xdc, impulse_response)
f2.xdc_excitation(xdc, excitation)
f2.xdc_focus_times(xdc, [0.0], np.mat([0.0,0.0]))

pos = np.array([0.0,0.0,15e-3]).reshape((1,3))
[hp, tstart0] = f2.calc_hp(xdc, pos)

fig  = plt.figure()

axes = [fig.add_subplot(121), fig.add_subplot(122)]
axes[0].plot(hp.T,'b')

[h, tstart1] = f2.calc_h(xdc, pos)
axes[1].plot(h.T,'b')

pos = pos.astype(np.float32)
a = fnm.ApertureFloat(nElements, width, kerf, height)

sp = a.sysparm
sp.timeDomainCalcType = fnm.TimeDomainCalcType.Plane
a.focus_type = fnm.FocusingType.Delays
a.sysparm = sp
a.fs = fs
a.normalize = True
a.c  = c
a.excitation       = [1.0]
a.impulse          = [1.0]
_tstart0, hp1 = a.CalcPwField(pos)
hp1 = np.roll(hp1,int((_tstart0 - tstart0)*fs-1.0))
axes[0].plot(hp1.T,'r')

a.excitation       = []
a.impulse          = []
_tstart1, h1 = a.CalcPwField(pos)
h1 = np.roll(h1,int((_tstart1 - tstart1)*fs-1.0))
axes[1].plot(h1.T,'r')
axes[0].legend(['SOFUS', 'Field II'])
axes[1].legend(['SOFUS', 'Field II'])
