#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 2; python-indent: 2; indent-tabs-mode: nil -*-

import sys
import argparse
import numpy as np

import unittest
import subprocess

import matplotlib.pyplot as plt

from timeit import default_timer as timer

import addpaths

import swig_fnm as fnm

from cw_acoustics import grating_lobe_angle, beam_pattern, accept_angle

def log_compress(pressure,dBrange=60):
  logp = np.abs(pressure)
  logp = logp / logp.max()
  logp = 20*np.log10(logp)
  logp[logp < -dBrange] = -dBrange
  logp[0,0] = 0
  logp[-1,-1] = -dBrange
  return logp

show = False
class FnmCwTest(unittest.TestCase):
  def setUp(self):
    """
    """
    f0 = 3.0e6              # excitation frequency [Hz]
    soundspeed = 1540       # [m/s]
    self.lamda = soundspeed / f0 # wavelength [m]

    #define a transducer structure/array
    nelex = 512
    neley = 1
    self.kerf = 5.0e-4

    self.width = 1.5*self.lamda # transducer width [m]
    height = 50e-3 # transducer height [m]

    ndiv = 40

    self.xdc = fnm.ApertureFloat(nelex, self.width, self.kerf, height)

    self.xdc.f0 = f0
    self.xdc.c  = soundspeed
    self.xdc.nDivW = ndiv
    self.xdc.nDivH = ndiv
    self.xdc.nthreads = 4

    # Setup grid
    pitch = self.width + self.kerf

    d = self.xdc.nelements * pitch

    xmin = -1.5 * d/2
    xmax = 1.5 * d/2
    ymin = 0
    ymax = 0
    zmin = 0.0
    zmax = 2*d

    nx = 170
    nz = 250

    dx = (xmax - xmin) / max(nx-1.0,1.0)
    dz = (zmax - zmin) / max(nz-1.0,1.0)
    self.xs0 = (np.r_[0:nx] - (nx-1.0)/2.0 +0.5) * dx
    self.zs0 = (np.r_[0:nz] + 0.5) * dz

    k = (2*np.pi)/self.lamda

    self.xs, self.zs = np.meshgrid(self.xs0, self.zs0, indexing='ij')

    self.nx = nx
    self.nz = nz

  def test_grating_lobe(self):
    pitch = self.width + self.kerf
    d = self.xdc.nelements * pitch

    apodization = self.xdc.apodization
    nActive = 120
    apodization[:] = 0.0
    apodization[self.xdc.nelements // 2 - nActive//2:self.xdc.nelements//2 + nActive//2] = 1.0
    self.xdc.apodization = apodization

    pos = np.c_[self.xs.flatten(), np.zeros(self.nx*self.nz), self.zs.flatten()].astype(np.float32)

    self.xdc.focus = [0,0,d]

    self.xdc.focus_type = fnm.FocusingType.Rayleigh

    start = timer()
    out = self.xdc.CalcCwFast(pos)[1]
    end = timer()
    print("Time elapsed: %f seconds" % (end-start))

    result3 = out.reshape((self.nx, self.nz))

    rho = 10000
    omega = 2*np.pi*self.xdc.f0
    pressure = -1j * omega * rho * np.exp(1j * omega * 0) * result3

    pressure = np.abs(pressure)
    pressure = pressure.astype(np.float64)
    logp = log_compress(pressure)

    steering = 0.0
    angles = grating_lobe_angle(steering, pitch, self.xdc.f0, c=1540)


    if show:
      plt.figure()
      plt.title('test: Grating Lobe')
      plt.imshow(logp,aspect='auto',extent=np.round(1000*np.r_[0,2*d,-1.5*d/2,1.5*d/2]),interpolation='none')
      plt.xlabel('Depth [mm]')
      plt.ylabel('Width [mm]')
      for angle in angles:
        dists = np.linspace(0, 500, 100)
        gx = dists*np.cos(angle)
        gy = dists*np.sin(angle)
        l, = plt.plot(gx, gy,'r-', linewidth=2.0)
      plt.legend([l], ['Grating lobe directions'])
    self.assertTrue(True)
  def test_acceptance_angle(self):
    from scipy.interpolate import RegularGridInterpolator

    # Symmetric
    nActive = 1
    pitch = self.width + self.kerf
    nElements = self.xdc.nelements
    d = nElements * pitch

    apodization = np.ones(self.xdc.nelements, dtype=np.float32)
    apodization[(nElements + nActive) // 2:] = 0
    apodization[:(nElements - nActive)//2] = 0

    self.xdc.apodization = apodization

    pos = np.c_[self.xs.flatten(),
                np.zeros(self.nx*self.nz),
                self.zs.flatten()].astype(np.float32)

    width = self.width

    for baffle in ['soft', 'hard']:
      sp = self.xdc.sysparm
      if baffle == 'soft':
        sp.soft_baffle = True
      else:
        sp.soft_baffle = False
      self.xdc.sysparm = sp

      out = self.xdc.CalcCwFast(pos)[1]
      result3 = out.reshape((self.nx, self.nz))

      rho = 10000
      omega = 2*np.pi*self.xdc.f0
      pressure = -1j * omega * rho * np.exp(1j * omega * 0) * result3

      pressure = np.abs(pressure)
      pressure = pressure.astype(np.float64)
      logp = log_compress(pressure)

      angle = accept_angle(0.5, width, self.xdc.f0, baffle_type=baffle)

      if show:
        plt.figure()
        plt.title('test: AcceptanceAngle')
        plt.imshow(logp, aspect='auto',
                   extent=np.round(1000*np.r_[0,2*d,-1.5*d/2,1.5*d/2]),
                   interpolation='none')
        plt.xlabel('Depth [mm]')
        plt.ylabel('Width [mm]')
        plt.colorbar()

      dists = np.linspace(0, 500, 100)
      gx = dists*np.cos(angle)
      dy = 3.0*d/2.0 / logp.shape[0]
      gy = dists*np.sin(angle) + (0.5 * pitch) / dy
      ls = []
      if show:
        ls.append(plt.plot(gx, gy,'r-', linewidth=2.0)[0])

      values = pressure
      f = RegularGridInterpolator((self.xs0, self.zs0), values)

      angles = np.linspace(-np.pi/4.0, np.pi/4.0, 100)

      bx = d*np.sin(angles)
      bz = d*np.cos(angles)
      bum = f(list(zip(bx,bz)))

      bum = bum / bum.max()

      iAngle = np.where(bum > 0.5)[0][-1]

      angle = angles[iAngle]
      gx = dists*np.cos(angle)
      dy = 3.0*d/2.0 / logp.shape[0]
      gy = dists*np.sin(angle) + (0.5 * pitch) / dy
      if show:
        ls.append(plt.plot(gx, gy,'r--', linewidth=2.0)[0])
        plt.legend(ls, [baffle + ' baffle (analytic)', 'simulation'])
    self.assertTrue(True)
  def test_beam_pattern(self):
    from scipy.interpolate import RegularGridInterpolator

    pitch = self.width + self.kerf
    nElements = self.xdc.nelements
    d = nElements * pitch

    # Symmetric
    nActive = nElements // 2
    apodization = np.ones(self.xdc.nelements,dtype=np.float32)
    apodization[(nElements + nActive) // 2:] = 0
    apodization[:(nElements - nActive)//2] = 0

    self.xdc.apodization = apodization

    self.xdc.focus = [0,0,d]

    pos = np.c_[self.xs.flatten(),
                np.zeros(self.nx*self.nz),
                self.zs.flatten()].astype(np.float32)

    out = self.xdc.CalcCwFast(pos)[1]

    result3 = out.reshape((self.nx, self.nz))

    nAngles = 100
    angles = np.linspace(-np.pi/4.0, np.pi/4.0, nAngles)

    bx = d*np.sin(angles)
    bz = d*np.cos(angles)

    values = np.real(result3)
    f = RegularGridInterpolator((self.xs0, self.zs0), values)
    rval = f(list(zip(bx,bz)))

    values = np.imag(result3)
    f = RegularGridInterpolator((self.xs0, self.zs0), values)
    ival = f(list(zip(bx,bz)))

    angPattern = beam_pattern(angles, pitch, nActive, self.xdc.f0,
                              kerf=self.kerf, weight=lambda x: np.cos(x)**2)
    angPattern = angPattern**2

    bum = rval + 1j*ival
    bum = bum * np.conj(bum)

    if show:
      plt.figure()
      plt.title('test: BeamPattern')
      l0, = plt.plot(angPattern/angPattern.max(),'r',linewidth=2.0)
      l1, = plt.plot(bum/bum.max(),'k')
      plt.legend([l0,l1], ['Analytic', 'FNM'])

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--show', type=str, choices=['False','True'], default='False', help='show the results')
  parser.add_argument('unittest_args', nargs='*')
  args = parser.parse_args()
  show = args.show == 'True'
  sys.argv[1:] = args.unittest_args
  unittest.main()

# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
