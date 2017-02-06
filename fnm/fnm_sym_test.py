#!/usr/bin/env python

import sys
import re
import numpy as np

import unittest
import subprocess

import matplotlib.pyplot as plt

import addpaths

import swig_fnm as fnm

from euler import euler2rot

plt.ion()

class SofusSymTest(unittest.TestCase):
  staticData = None
  def setUp(self):
    width  = 1.0
    kerf   = 0.0
    height = 2.0
    
    self.a = fnm.ApertureFloat(1, width, kerf, height)
    self.a.c          = 1.0
    self.a.f0         = 1.0

    pos = np.r_[3,2,5]
    self.pos = pos.reshape([1,3]).astype(np.float32)

    self.a.focus = self.pos.flatten()
    
    retval, self.ref = self.a.CalcCwField(self.pos)

  def rotate_axes(self, limits):
    retval = True
    maxdiff = 0.0
    for alpha in np.r_[limits[0,0]:limits[0,1]:limits[0,2]]:
      for beta in np.r_[limits[1,0]:limits[1,1]:limits[1,2]]:
        for gamma in np.r_[limits[2,0]:limits[2,1]:limits[2,2]]:
          rmat = euler2rot(alpha,beta,gamma,conv='yxy')
          pos2 = np.dot(rmat,self.pos.T).T
          pos2 = pos2.astype(np.float32)
          
          elements = self.a.elements
          elements[0,5:8] = [alpha,beta,gamma]
          self.a.elements = elements
    
          self.a.focus = pos2.flatten()
          _, hp2 = self.a.CalcCwField(pos2)
          hp2 = hp2.flatten()
          diff = np.abs(self.ref-hp2) / (np.abs(self.ref)+np.abs(hp2))
          maxdiff = max(maxdiff, diff[0])
          retval = retval and (diff[0] < 0.05)
    if not(retval):
      print('max. diff: %f' % maxdiff)
    return retval

  def test_triple(self):
    limits = np.c_[[0.0,       0.0,       0.0       ],
                   [np.pi/2.0, np.pi/2.0, np.pi/2.0 ],
                   [np.pi/10.0,np.pi/10.0,np.pi/10.0]]

    # Fails (when gamma is included)
    success = self.rotate_axes(limits)
    self.assertTrue(success)
  
  def test_double(self):
    limits = np.c_[[0.0,       0.0,       0.0],
                   [np.pi/2.0, np.pi/2.0, 0.2],
                   [np.pi/20.0,np.pi/20.0,0.2]]

    success = self.rotate_axes(limits)
    self.assertTrue(success)

    success = self.rotate_axes(limits[[2,0,1]])
    self.assertTrue(success)

    # Fails (when gamma is included)
    success = self.rotate_axes(limits[[0,2,1]])
    self.assertTrue(success)
    
  def test_single(self):
    limits = np.c_[[0.0,       0.0,   0.0],
                   [np.pi/2.0, 0.2,   0.2],
                   [np.pi/20.0,0.2,   0.2]]

    success = self.rotate_axes(limits)
    self.assertTrue(success)

    limits = np.c_[[0.0,   0.0,       0.0],
                   [0.2,   np.pi/2.0, 0.2],
                   [0.2,   np.pi/20.0,0.2]]

    success = self.rotate_axes(limits)
    self.assertTrue(success)

    limits = np.c_[[0.0,   0.0,       0.0       ],
                   [0.2,   0.2,       np.pi/2.0 ],
                   [0.2,   0.2,       np.pi/20.0]]

    success = self.rotate_axes(limits)
    self.assertTrue(success)
    
if __name__ == '__main__':
  unittest.main()
    
# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# py-indent-offset: 2 #
# indent-tabs-mode: nil #
# End: #
    
