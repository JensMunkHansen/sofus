#!/usr/bin/env python

import os, sys
import argparse
import re
import numpy as np

import unittest
import subprocess

import types

import matplotlib.pyplot as plt

import addpaths

from euler import euler2rot

import swig_fnm as fnm

from reference import rect
from reference import pulse


plt.ion()

show = False
class FnmTransientTest(unittest.TestCase):
  def setUp(self):
    """
    """
    width = height = 1.0e-3
    self.a = fnm.ApertureFloat(1,width,0.0,height)

  #@unittest.skip("Crashes")
  def test_something(self):
    width = height = 5.0e-3
    depth = 20e-3
    f0 = 1e6
    c = 1500.0
    ndiv = 8
    fs = 10e6

    a = fnm.ApertureFloat(1,width,0.0,height)

    a.f0 = f0
    a.c = c
    a.nDivW = ndiv
    a.nDivH = ndiv

    # Calculate the pressure
    a.fs = fs # ISSUE

    nCycles = 3.0

    a.w = nCycles / a.f0
    a.focus_type = fnm.FocusingType.Pythagorean

    limits = np.c_[[0.0,                   0.0,       0.0       ],
                   [np.pi/2.0,       np.pi/2.0,       2*np.pi   ],
                   [np.pi/2.0,       np.pi/2.0,       np.pi/2.0 ]]

    # Lower for edge 0
    #pos = np.r_[np.c_[depth, depth, depth]].astype(np.float32)
    # Upper for edge 0
    pos = np.r_[np.c_[-depth, 0.0, depth]].astype(np.float32)

    if show:
      plt.figure()

    #res = a.CalcTransientSingleElementNoDelay(pos, 0x02)

    if 0:
      for alpha in np.r_[limits[0,0]:limits[0,1]:limits[0,2]]:
        for beta in np.r_[limits[1,0]:limits[1,1]:limits[1,2]]:
          for gamma in np.r_[limits[2,0]:limits[2,1]:limits[2,2]]:
            rmat = euler2rot(alpha,beta,gamma,conv='yxz')
            pos2 = np.dot(rmat,pos.T).T
            pos2 = pos2.astype(np.float32)
            res = a.CalcTransientSingleElementNoDelay(pos2)
            if show:
              plt.plot(res[1].T)

  #@unittest.skip("")
  def test_against_reference(self):
    hh=7.5e-3
    hw=7.5e-3
    fs = 20e6
    f0 = 1e6
    c = 1500
    nCycles = 3
    ndiv = 4

    r = rect(hh=hh,hw=hw, ndiv=ndiv)
    p = pulse(fs=fs,f0=f0,c=c,nCycles=nCycles)

    points = np.r_[np.c_[0.0, 8.3e-3, 8e-3],
                   np.c_[8.3e-3, 0.0, 8e-3]]

    xdc = fnm.ApertureFloat(1,2*hw,0.0,2*hh)
    xdc.fs = fs
    xdc.f0 = f0
    xdc.fc = f0
    xdc.excitation_type = fnm.ExcitationType.ExcitationTypeToneBurst
    xdc.w = nCycles / xdc.f0
    xdc.c = c
    xdc.nDivW = ndiv
    xdc.nDivH = ndiv
    xdc.focus_type = fnm.FocusingType.Pythagorean

    if show:
      fh, axes = plt.subplots(2,3)
      axes = axes.flatten()

    for j in range(2):
      point = points[j:j+1,:]
      ref = r.calc_direct_response(p,point)
      t0, test = xdc.CalcPwFnmThreaded(point.astype(np.float32), 0x01)
      test = test.flatten()
      test = np.pad(test,(0,len(ref)),mode='constant')
      test = np.roll(test, int(t0*fs))
      if show:
        axes[0].plot(test)
        axes[0].plot(ref)

      for i in range(1,5):
        ref = r.calc_edge_response(p,i-1,point)
        t0, test = xdc.CalcPwFnmThreaded(point.astype(np.float32), 1<<i)
        test = test.flatten()
        test = np.pad(test, (0,len(ref)), mode='constant')
        test = np.roll(test, int(t0*fs))
        if show:
          axes[i].plot(test)
          axes[i].plot(ref)


  #@unittest.skip("")
  def test_symmetry(self):
    """
    Test symmetry when rotating both point and element
    """
    width = height = 1.0e-3
    depth = 1e-3
    f0 = 1e6
    c = 1500.0
    ndiv = 4
    fs = 4 * f0

    a = fnm.ApertureFloat(1,width,0.0,height)

    a.f0 = f0
    a.c = c
    a.nDivW = ndiv
    a.nDivH = ndiv

    # Calculate the pressure
    a.fs = fs

    nCycles = 3.0

    a.w = nCycles / a.f0

    limits = np.c_[[0.0,       0.0,       0.0       ],
                   [np.pi/2.0,       np.pi/2.0,       2*np.pi   ],
                   [np.pi/2.0,       np.pi/2.0,       np.pi/2.0 ]]

    pos = np.r_[np.c_[0.0, 1.5*width, depth]].astype(np.float32)

    if show:
      plt.figure()
    for alpha in np.r_[limits[0,0]:limits[0,1]:limits[0,2]]:
      for beta in np.r_[limits[1,0]:limits[1,1]:limits[1,2]]:
        for gamma in np.r_[limits[2,0]:limits[2,1]:limits[2,2]]:
          rmat = euler2rot(alpha,beta,gamma,conv='yxz')
          pos2 = np.dot(rmat,pos.T).T
          pos2 = pos2.astype(np.float32)
          res = a.CalcTransientSingleElementNoDelay(pos2)
          if show:
            plt.plot(res[1].T)

    pos = np.r_[np.c_[1.5*width, 1.5*width, depth]].astype(np.float32)

    if show:
      plt.figure()
    for alpha in np.r_[limits[0,0]:limits[0,1]:limits[0,2]]:
      for beta in np.r_[limits[1,0]:limits[1,1]:limits[1,2]]:
        for gamma in np.r_[limits[2,0]:limits[2,1]:limits[2,2]]:
          rmat = euler2rot(alpha,beta,gamma,conv='yxz')
          pos2 = np.dot(rmat,pos.T).T
          pos2 = pos2.astype(np.float32)
          res = a.CalcTransientSingleElementNoDelay(pos2)
          if show:
            plt.plot(res[1].T)




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
