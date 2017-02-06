#!/usr/bin/env python
import sys
import re
import numpy as np

import unittest
import subprocess

import addpaths

import swig_fnm as fnm
from dicts import dotdict

class FnmReferenceTest(unittest.TestCase):
    def setUp(self):
      # Tolerance between SIMD and C reference implementation
      self.eps = 1000.0 * np.finfo(np.float32).eps

    def test_inside(self):
      # SSE with C refernece
      retval = self.compareWithReference(location='inside', method1='CalcCwFast', eps=self.eps)
      self.assertEqual(retval[0],True)

    def test_outside(self):
      # SSE with C refernece
      retval = self.compareWithReference(location='outside',method1='CalcCwFast',eps=self.eps)
      self.assertEqual(retval[0],True)

    @staticmethod
    def compareWithReference(**kwargs):
      opt = dotdict({'method0'  : 'CalcCwFieldRef',
                     'method1'  : 'CalcCwFieldRef',
                     'location' : 'inside',
                     'ndiv'     : 2,
                     'eps'      : 1e-5,
                     'quiet'    : True})

      opt.update(**kwargs)
      quiet = opt.quiet
      
      areas = [2.0,3.0,4.0,5.0]
      widths  = np.array([areas[0],1.0,areas[2],1.0],dtype=np.float32)
      heights = np.array([1.0,areas[1],1.0,areas[3]],dtype=np.float32)
      
      xsign = np.array([1.0,-1.0,-1.0,1.0],dtype=np.float32)
      ysign = np.array([1.0,1.0,-1.0,-1.0],dtype=np.float32)
      
      # Ensure that we either slightly inside or outside
      scale = 40.0 * np.finfo(np.float32).eps
      
      epss = dict({'inside'  : -scale,
                   'outside' :  scale})
      
      ndiv = opt.ndiv
      iCorners = [0,1,2,3]
      
      eps    = epss[opt.location]
      method0 = opt.method0
      method1 = opt.method1
      
      scatters = np.c_[(widths+eps)/2.0 * xsign,
                       (heights+eps)/2.0 * ysign,
                       np.ones(4,dtype=np.float32)]
      if not(quiet):
        print('Testing: %s vs %s: %s' % (opt.method0,opt.method1, opt.location))
      success = True
      for iCorner in iCorners:
        a = fnm.ApertureFloat(1,
                              float(widths[iCorner]),
                              float(0.0),
                              float(heights[iCorner]))
        a.nthreads = 1
        a.nDivH = ndiv
        a.nDivW = ndiv
        a.f0 = 1
        a.c  = 1
        
        exec('result0 = a.'+method0+'([scatters[iCorner]])[1].flatten()')
        exec('result1 = a.'+method1+'([scatters[iCorner]])[1].flatten()')
        diff = np.abs(result0)-np.abs(result1)
        reldiff = diff / np.mean(np.abs(result0)+np.abs(result1))
        success = success and (reldiff < opt.eps)
        if not(quiet):
          print('Relative difference: %f' % (reldiff))
      return success

if __name__ == '__main__':
    unittest.main()

