#!/usr/bin/env python
import sys
import re
import argparse
import numpy as np

import unittest
import subprocess

import addpaths

from fnm import rect
import swig_fnm as fnm
from dicts import dotdict

import matplotlib.pyplot as plt


show = False
quiet = True

plt.ion()

class FnmReferenceTest(unittest.TestCase):
  def setUp(self):
    # Tolerance between Python and C++ implementation
    self.eps = 3.0*np.finfo(np.float32).eps # Linux works 1.0
  def test_inside(self):
    retval = self.compareWithPython(location='inside', method='CalcCwFieldRef', eps=self.eps)
    self.assertEqual(retval[0],True)
  def test_outside(self):
    retval = self.compareWithPython(location='outside',method='CalcCwFieldRef', eps=self.eps)
    self.assertEqual(retval[0],True)

  @staticmethod
  def compareWithPython(**kwargs):
    opt = dotdict({'method'   : 'CalcCwFieldRef',
                   'location' : 'inside',
                   'eps'      : 1e-5})
    opt.update(**kwargs)
    
    areas = [2.0,3.0,4.0,5.0]
    widths  = np.array([areas[0],1.0,areas[2],1.0],dtype=np.float32)
    heights = np.array([1.0,areas[1],1.0,areas[3]],dtype=np.float32)
    
    xsign = np.array([1.0,-1.0,-1.0,1.0],dtype=np.float32)
    ysign = np.array([1.0,1.0,-1.0,-1.0],dtype=np.float32)
    
    # Ensure that we either slightly inside or outside
    scale = 0.1
    
    offsets = dict({'inside'  : -scale,
                    'outside' :  scale})
    
    ndiv = 2
    iCorners = [0,1,2,3]
    
    offset    = offsets[opt.location]
    method = opt.method
    
    scatters = np.c_[(widths+offset)/2.0 * xsign,
                     (heights+offset)/2.0 * ysign,
                     np.ones(4,dtype=np.float32)]

    if not(quiet):
      print('Testing: %s vs %s: projection %s element' % (opt.method,'H_accurate (Python)',opt.location))

    success = True
    for iCorner in iCorners:
      a = fnm.ApertureFloat(1,float(widths[iCorner]),
                            float(0.0),
                            float(heights[iCorner]))
      a.nthreads = 1
      a.nDivH = ndiv
      a.nDivW = ndiv
      a.f0 = 1
      a.c  = 1
      if not(quiet):
        if show:
          a.show()
          ax = plt.gca()
          ax.scatter([scatters[iCorner][0]],
                     [scatters[iCorner][1]],
                     [scatters[iCorner][2]],color='black',marker='o',s=20)
          
          ax.set_xlim([-widths.max(),widths.max()])
          ax.set_ylim([-heights.max(),heights.max()])
          ax.set_zlim([-5,5])
          ax.set_aspect('equal')
          plt.show()
    
      exec('result = a.'+method+'([scatters[iCorner]])[1].flatten()')
    
      r = rect(hw=widths[iCorner]/2.0,hh=heights[iCorner]/2.0,center=[0,0,0],nAbcissa=ndiv)
      
      k = 2*np.pi / (a.c / a.f0)
    
      xs = np.ones((1,1))*scatters[iCorner][0]
      ys = np.ones((1,1))*scatters[iCorner][1]
      zs = np.ones((1,1))*scatters[iCorner][2]
    
      reference = r.H_accurate(xs,ys,zs,k).flatten()
      assert(np.abs(reference) - np.abs(result) < np.finfo(np.float32).eps)
      diff = np.abs(reference)-np.abs(result)
      reldiff = diff / np.mean(np.abs(reference)+np.abs(result))
      nextSuccess = reldiff < opt.eps
      if not(quiet) or not(nextSuccess):
        print('Relative difference: %g' % (reldiff))
      success = success and nextSuccess
    return success

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--show', type=str, choices=['False','True'], default='False', help='show the results')
  parser.add_argument('--quiet', type=str, choices=['False','True'], default='True', help='show the results')
  parser.add_argument('unittest_args', nargs='*')
  args = parser.parse_args()
  show = args.show == 'True'
  quiet = args.quiet == 'True'
  sys.argv[1:] = args.unittest_args
  unittest.main()

# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# py-indent-offset: 2 #
# indent-tabs-mode: nil #
# End: #
